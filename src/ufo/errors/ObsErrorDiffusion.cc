
/*
 * (C) Copyright 2025 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <cmath>
#include <numeric>
#include <unordered_map>
#include <utility>
#include <vector>

#include "atlas/functionspace/NodeColumns.h"
#include "atlas/grid/Grid.h"
#include "atlas/mesh/actions/BuildCellCentres.h"
#include "atlas/mesh/actions/BuildHalo.h"
#include "atlas/mesh/actions/BuildXYZField.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/MeshBuilder.h"
#include "atlas/meshgenerator/MeshGenerator.h"
#include "atlas/output/Gmsh.h"
#include "atlas/util/Geometry.h"
#include "atlas/util/KDTree.h"

#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"

#include "oops/assimilation/GMRESR.h"
#include "oops/base/FieldSet3D.h"
#include "oops/base/IdentityMatrix.h"
#include "oops/mpi/Scope.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FunctionSpaceHelpers.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/errors/ObsErrorDiffusion.h"

namespace ufo {

// -----------------------------------------------------------------------------

static ObsErrorMaker<ObsErrorDiffusion> makerDiffusion_("diffusion");

// -----------------------------------------------------------------------------

ObsErrorDiffusion::ObsErrorDiffusion(const Parameters_ & params,
                                     ioda::ObsSpace & odb,
                                     const eckit::mpi::Comm &timeComm)
  : ObsErrorBase(timeComm),
    params_(params),
    stddev_(odb, "ObsError"),
    inverseVariance_(odb),
    comm_(odb.comm())
{
  oops::Log::trace() << "ObsErrorDiffusion::ObsErrorDiffusion constructor start" << std::endl;
  ASSERT(params_.var.value().size() == 1);  // only one variable for now
  inverseVariance_ = stddev_;
  inverseVariance_ *= stddev_;
  inverseVariance_.invert();
  oops::Log::trace() << "ObsErrorDiffusion::ObsErrorDiffusion constructed" << std::endl;
}

// -----------------------------------------------------------------------------

// Refresh the obs-error standard deviations to \p obserr and, on the first call, build
// everything the correlation operator C needs: the global diffusion mesh (from all obs
// across all ranks), the oops::Diffusion operator, the obs->mesh-node map, and the
// normalization factors Gamma. The mesh and normalization depend only on the obs
// locations, QC mask, and length scale -- not on the error values -- so later calls just
// refresh the standard deviations and reuse the cached build. R is subsequently applied
// by multiply() / inverseMultiply() / randomize().
void ObsErrorDiffusion::update(const ioda::ObsVector & obserr) {
  oops::Log::trace() << "ObsErrorDiffusion update: start " << std::endl;

  // ------------------------------------------------------------
  // Refresh the error standard deviations and inverse variance
  // ------------------------------------------------------------
  stddev_ = obserr;
  inverseVariance_ = stddev_;
  inverseVariance_ *= stddev_;
  inverseVariance_.invert();

  // ------------------------------------------------------------
  // Collect this rank's valid (QC-passing) obs
  // ------------------------------------------------------------
  // Obs are keyed on their source-file index (ObsSpace::index()), which is
  // partition-invariant, so every rank places a given obs at the same mesh node and R
  // is the same operator for any PE count; globalUniqueConsecutiveLocationIndex is not
  // (it renumbers by rank, perturbing the stochastic normalization). Requires a
  // non-overlapping, geographically contiguous distribution (e.g. Atlas)

  // *for single channel/var* all obs locations on this PE (including obs not passing qc)
  const int nlocs = obserr.nlocs();
  // *for single channel/var* nobs is number of obs passing QC across all mpi tasks
  const int nobs = obserr.nobs();
  const int nparts = comm_.size();
  std::vector<float> lons(nlocs);
  std::vector<float> lats(nlocs);
  obserr.space().get_db("MetaData", "longitude", lons);
  obserr.space().get_db("MetaData", "latitude", lats);

  const double missing = util::missingValue<double>();

  const ioda::Distribution & dist = *obserr.space().distribution();
  if (nparts > 1 && !dist.isNonoverlapping()) {
    throw eckit::BadParameter("ObsErrorDiffusion requires a non-overlapping distribution "
        "(e.g. Atlas) on more than one PE; '" + dist.name() + "' is overlapping.");
  }
  if (nparts > 1 && dist.name() != "Atlas" && !params_.allowAnyDistribution.value()) {
    throw eckit::BadParameter("ObsErrorDiffusion requires a geographically contiguous "
        "distribution on more than one PE, and '" + dist.name() + "' is not guaranteed to be "
        "one. 'Atlas' is the only distribution that guarantees it; others may be contiguous "
        "depending on the setup, in which case set 'allow any distribution: true'.");
  }

  const std::vector<std::size_t> & srcIndex = obserr.space().index();
  std::vector<std::size_t> myGidx;     // source-file location index per valid obs
  std::vector<double> myLon, myLat;
  std::vector<int> myLocalIdx;         // local ObsVector index per valid obs
  myGidx.reserve(nlocs);
  myLon.reserve(nlocs);
  myLat.reserve(nlocs);
  myLocalIdx.reserve(nlocs);
  for (int i = 0; i < nlocs; ++i) {
    if (obserr[i] == missing) continue;  // passivated/QC-failed obs are set to missingValue
    myGidx.push_back(srcIndex[i]);
    myLon.push_back(lons[i]);
    myLat.push_back(lats[i]);
    myLocalIdx.push_back(i);
  }

  // ------------------------------------------------------------
  // If update() has been already called before:
  // ------------------------------------------------------------
  // TEMPORARY: the reuse guard only checks this rank's obs count is unchanged. That
  // holds within one minimization but not across outer loops whose QC swaps obs at
  // equal count (would reuse a stale mesh). Doing that right means deciding how to
  // handle "missing" obs between meshes; revisit when multi-outer-loop QC changes are
  // needed.
  if (meshAndNormBuilt_) {
    ASSERT(static_cast<int>(myGidx.size()) == builtObsCount_);
    oops::Log::info() << "ObsErrorDiffusion::update: reusing diffusion mesh and "
        "normalization (obs set unchanged)" << std::endl;
    oops::Log::trace() << "ObsErrorDiffusion update: end (reused)" << std::endl;
    return;
  }

  // ------------------------------------------------------------
  // Gather all valid obs across ranks, in a canonical order
  // ------------------------------------------------------------
  // allGatherv concatenates in rank order; the owner of each gathered entry is
  // recovered from the per-rank counts. The rank-ordered intermediates are scoped so
  // these global-sized arrays are freed before the (memory-hungry) mesh build.
  std::vector<atlas::PointLonLat> obsnodes;  // canonical-order obs coords
  std::vector<int> obsOwner;                 // canonical-order owning rank
  std::vector<std::size_t> sortedGidx;       // canonical order = ascending source-file index
  {
    std::vector<int> counts(nparts, 0);
    const int myCount = static_cast<int>(myGidx.size());
    comm_.allGather(myCount, counts.begin(), counts.end());

    std::vector<std::size_t> gGidx = myGidx;  oops::mpi::allGatherv(comm_, gGidx);
    std::vector<double>      gLon  = myLon;   oops::mpi::allGatherv(comm_, gLon);
    std::vector<double>      gLat  = myLat;   oops::mpi::allGatherv(comm_, gLat);
    ASSERT(static_cast<int>(gGidx.size()) == nobs);

    std::vector<int> gOwner(gGidx.size());
    for (int r = 0, pos = 0; r < nparts; ++r)
      for (int k = 0; k < counts[r]; ++k) gOwner[pos++] = r;

    // Canonical, partition-invariant ordering: sort by source-file index.
    std::vector<std::size_t> order(gGidx.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
              [&](std::size_t a, std::size_t b) { return gGidx[a] < gGidx[b]; });

    obsnodes.resize(gGidx.size());
    obsOwner.resize(gGidx.size());
    sortedGidx.resize(gGidx.size());
    for (std::size_t k = 0; k < gGidx.size(); ++k) {
      const std::size_t s = order[k];
      atlas::PointLonLat p(gLon[s], gLat[s]);
      p.normalise();  // to [0, 360), the convention used by the mesh and control grid
      obsnodes[k] = p;
      obsOwner[k] = gOwner[s];
      sortedGidx[k] = gGidx[s];
    }
  }
  const std::size_t nGlobalObs = obsnodes.size();

  oops::Log::trace() << "Constructing diffusion mesh with " << nobs
                     << " valid observations across " << nparts << " PE(s)" << std::endl;

  // ------------------------------------------------------------
  // Optional coarse control grid
  // ------------------------------------------------------------
  // Built from the GLOBAL obs set so it is identical on every rank. Control points
  // carry no obs data; each inherits the owning rank of its nearest obs so the
  // partition stays geographically consistent. Not enabled by default -- grid spacing
  // and remove-within must be set explicitly for the observation network geometry.
  std::vector<atlas::PointLonLat> ctrlPoints;
  std::vector<int> ctrlOwner;
  if (params_.controlGrid.value()) {
    const int gridSpacing = params_.controlGrid.value()->gridSpacing.value();
    const double removeWithin = params_.controlGrid.value()->removeWithin.value();
    createControlGrid(obsnodes, obsOwner, removeWithin, gridSpacing, ctrlPoints, ctrlOwner);
  }
  // obsOffset_ = number of control points; control points occupy global grid
  // indices [0, obsOffset_), obs occupy [obsOffset_, obsOffset_ + nGlobalObs).
  obsOffset_ = static_cast<atlas::idx_t>(ctrlPoints.size());

  // ------------------------------------------------------------
  // Build the distributed Delaunay mesh and diffusion operator
  // ------------------------------------------------------------
  // The triangulation is over the full global point set (control points + obs, built
  // identically on every rank); the Distribution carries IODA's per-obs ownership so
  // each rank extracts its owned nodes plus a 1-ring ghost layer. On 1 PE this is a
  // single partition with no halo (the serial case).
  // NOTE: holding and triangulating the global point set on every rank is O(nobs)
  // memory/work per rank, which won't scale to very large obs counts -- a known
  // limitation, good enough for now. It is mostly transient: the global-sized
  // intermediates are freed at last use (see the scoping and early frees), and after
  // update() only this rank's mesh partition (+ halo) remains -- plus the global point
  // coordinates, which atlas keeps attached to the mesh (Mesh::grid()).
  const int npGrid = static_cast<int>(ctrlPoints.size() + nGlobalObs);
  std::vector<atlas::PointXY> pts;
  std::vector<int> ownerArr;
  pts.reserve(npGrid);
  ownerArr.reserve(npGrid);
  for (std::size_t k = 0; k < ctrlPoints.size(); ++k) {
    pts.emplace_back(ctrlPoints[k].lon(), ctrlPoints[k].lat());
    ownerArr.push_back(ctrlOwner[k]);
  }
  for (std::size_t k = 0; k < nGlobalObs; ++k) {
    pts.emplace_back(obsnodes[k].lon(), obsnodes[k].lat());
    ownerArr.push_back(obsOwner[k]);
  }
  // global-sized and no longer needed; free before the mesh build
  std::vector<atlas::PointLonLat>().swap(obsnodes);
  std::vector<int>().swap(obsOwner);

  // The generator, Distribution and NodeColumns use eckit's default communicator;
  // scope it to comm_ for the build (oops::mpi::Scope restores the previous default
  // on exit, including on exceptions).
  atlas::Mesh mesh;
  atlas::functionspace::NodeColumns fspace;
  {
    const oops::mpi::Scope commScope(comm_);
    atlas::UnstructuredGrid grid(std::move(pts));
    atlas::grid::Distribution distribution(nparts, npGrid, ownerArr.data());
    mesh = atlas::MeshGenerator("delaunay").generate(grid, distribution);
    if (nparts > 1) {
      // 1-ring ghost layer for the nearest-neighbour stencil
      atlas::mesh::actions::build_halo(mesh, 1);
    }
    // oops::Diffusion::calculateDerivedGeom requires the node "xyz" and cell "centre"
    // fields. The delaunay generator builds them on the global mesh, but does not
    // carry them onto the per-rank mesh, so build them here (xyz first;
    // BuildCellCentres needs it).
    atlas::mesh::actions::BuildXYZField()(mesh);
    atlas::mesh::actions::BuildCellCentres()(mesh);
    fspace = atlas::functionspace::NodeColumns(mesh);
  }

  // Owned-points mask (inverse of atlas ghost), as oops::GeometryData expects.
  atlas::FieldSet fset;
  atlas::Field owned = fspace.createField<int>(atlas::option::name("owned")
                                               | atlas::option::levels(1));
  {
    auto ownedView = atlas::array::make_view<int, 2>(owned);
    auto ghostView = atlas::array::make_view<int, 1>(fspace.ghost());
    for (atlas::idx_t i = 0; i < ghostView.shape(0); ++i)
      ownedView(i, 0) = (ghostView(i) > 0 ? 0 : 1);
  }
  fset.add(owned);

  geom_.reset(new oops::GeometryData(fspace, fset, true, comm_));
  diffusionGeom_ = oops::Diffusion::calculateDerivedGeom(std::as_const(*geom_));
  diffusion_.reset(new oops::Diffusion(std::as_const(*geom_), diffusionGeom_));

  const atlas::idx_t npts = geom_->functionSpace().size();

  // ------------------------------------------------------------
  // Map this rank's valid obs -> local mesh node
  // ------------------------------------------------------------
  // atlas assigns mesh node global_index = (global grid index) + 1; obs at global
  // grid index (obsOffset_ + canonical position) therefore has global_index
  // (obsOffset_ + pos + 1). Owned obs are guaranteed present (non-ghost) on this rank.
  auto meshGidxView = atlas::array::make_view<atlas::gidx_t, 1>(mesh.nodes().global_index());
  std::unordered_map<atlas::gidx_t, atlas::idx_t> gnode2node;
  gnode2node.reserve(npts);
  for (atlas::idx_t i = 0; i < npts; ++i) gnode2node[meshGidxView(i)] = i;

  // localObs2node_ is indexed by LOCAL ObsVector location: the owned mesh node for
  // each valid obs, and -1 for QC-failed locations, which multiply()/randomize() skip.
  localObs2node_.assign(nlocs, -1);
  for (std::size_t kk = 0; kk < myGidx.size(); ++kk) {
    // canonical position = rank of this obs's source-file index in the sorted global set
    const auto gIt = std::lower_bound(sortedGidx.begin(), sortedGidx.end(), myGidx[kk]);
    ASSERT(gIt != sortedGidx.end() && *gIt == myGidx[kk]);
    const std::size_t pos = gIt - sortedGidx.begin();
    const atlas::gidx_t meshG = static_cast<atlas::gidx_t>(obsOffset_ + pos) + 1;
    const auto it = gnode2node.find(meshG);
    ASSERT(it != gnode2node.end());
    localObs2node_[myLocalIdx[kk]] = it->second;
  }
  // global-sized and no longer needed; free before the normalization iterations
  std::vector<std::size_t>().swap(sortedGidx);

  // ------------------------------------------------------------
  // Optional: save mesh + point-type marker for gmsh visualization
  // ------------------------------------------------------------
  if (params_.outputDiffusionMesh.value()) {
    // Create marker field: 0 = control grid point, 1 = obs point
    atlas::Field markerField = geom_->functionSpace().createField<double>(
        atlas::option::levels(1) | atlas::option::name("pointType"));
    auto v_marker = atlas::array::make_view<double, 2>(markerField);
    v_marker.assign(0.0);

    // Obs nodes are global grid indices >= obsOffset_, i.e. global_index > obsOffset_.
    for (atlas::idx_t i = 0; i < npts; ++i) {
        if (meshGidxView(i) > obsOffset_) v_marker(i, 0) = 1.0;
    }
    fset.add(markerField);  // Add to fieldset before writing
    // Write mesh with the marker field
    const std::string filename = "diffusion_mesh.msh";
    atlas::output::Gmsh gmsh(filename,
        atlas::util::Config("coordinates", "xyz")
        | atlas::util::Config("ghost", true));  // enables viewing halos per task
    gmsh.write(mesh);
    gmsh.write(fset, geom_->functionSpace());
    oops::Log::info() << "ObsErrorDiffusion: Saved diffusion mesh to "
                        << filename << std::endl;
  }

  // ------------------------------------------------------------
  // Create horizontal length scales for normalization
  // ------------------------------------------------------------
  atlas::Field hzScales = geom_->functionSpace().createField<double>(
      atlas::option::levels(1) | atlas::option::name("hzScales"));
  auto v_hzScales = atlas::array::make_view<double, 2>(hzScales);
  // `correlation lengthscale` is a Gaspari-Cohn cutoff radius by default; the
  // 1/3.67 factor converts it to the Daley length oops::Diffusion expects.
  // `as gaussian: true` instead passes the value through as a Daley length.
  const double rawLscale = params_.lscale;
  const bool asGaussian = params_.asGaussian;
  const double daleyLscale = asGaussian ? rawLscale : (rawLscale / 3.67);
  v_hzScales.assign(daleyLscale);
  // Set the scales for control points to 0 so they don't influence the diffusion of obs points.
  // Control points are global grid indices [0, obsOffset_), i.e. global_index in [1, obsOffset_].
  for (atlas::idx_t i = 0; i < npts; ++i) {
    if (meshGidxView(i) <= obsOffset_) v_hzScales(i, 0) = 0;
  }
  atlas::FieldSet scales;
  scales.add(hzScales);

  oops::Log::info() << "  correlation lengthscale: " << rawLscale << " m ("
                    << (asGaussian ? "Daley length" : "Gaspari-Cohn cutoff")
                    << ") -> Daley length " << daleyLscale << " m" << std::endl;
  diffusion_->setParameters(scales);

  // ------------------------------------------------------------
  // Calculate horizontal normalization coefficients (Γ)
  // ------------------------------------------------------------
  const size_t levels = 1;
  const int randomizationIterations = params_.normalizationIterations;

  oops::Log::info() << "Calculating horizontal normalization using "
                  << params_.normalizationIterations << " iterations" << std::endl;

  // Field to hold normalization factors Γ
  this->hzNorm_ = geom_->functionSpace().createField<double>(
      atlas::option::levels(levels) | atlas::option::name("hzNorm"));
  auto v_norm = atlas::array::make_view<double, 2>(hzNorm_);
  v_norm.assign(1.0);

  // Temporary fields for running mean (m) and variance sum (s)
  atlas::Field m = geom_->functionSpace().createField<double>(atlas::option::levels(levels));
  atlas::Field s = geom_->functionSpace().createField<double>(atlas::option::levels(levels));

  auto v_m = atlas::array::make_view<double, 2>(m);
  auto v_s = atlas::array::make_view<double, 2>(s);

  v_m.assign(0.0);
  v_s.assign(0.0);

  // Main randomization loop
  for (int itr = 1; itr <= randomizationIterations; ++itr) {
    if (itr % std::max(1, randomizationIterations/10) == 0) {
      oops::Log::info() << "  progress: " << (100*itr/randomizationIterations)
                        << "%" << std::endl;
    }

    // ------------ Generate random field ------------
    // only one field needed since all variables on
    // same diffusion mesh will use same normalization
    atlas::FieldSet randSet = util::createRandomFieldSet(
      geom_->comm(), geom_->functionSpace(),
      std::vector<size_t>{levels}, std::vector<std::string>{"rand"});

    // ---- Apply sqrt diffusion ----
    diffusion_->multiplySqrtTL(randSet, oops::Diffusion::Mode::HorizontalOnly);

    // ---- Update Welford running variance ----
    auto v_rand = atlas::array::make_view<double, 2>(randSet["rand"]);

    for (atlas::idx_t i = 0; i < npts; i++) {
      double f = v_rand(i, 0);
      double old_m = v_m(i, 0);
      double new_m = old_m + (f - old_m) / itr;
      v_s(i, 0) += (f - old_m) * (f - new_m);
      v_m(i, 0) = new_m;
    }
  }

  // ---- Finalize normalization: Γ_i = 1/sqrt(Var_i) ----
  for (atlas::idx_t i = 0; i < npts; i++) {
    if (v_s(i, 0) > 0.0) {
      v_norm(i, 0) = 1.0 / std::sqrt(v_s(i, 0) / (randomizationIterations - 1));
    }
  }
  geom_->functionSpace().haloExchange(hzNorm_);

  // ------------------------------------------------------------
  // Record the build so later update() calls hit the reuse guard
  // ------------------------------------------------------------
  builtObsCount_ = static_cast<int>(myGidx.size());
  meshAndNormBuilt_ = true;

  oops::Log::trace() << "ObsErrorDiffusion update: end" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsErrorDiffusion::multiply(ioda::ObsVector & dy) const {
  // R * dy = D^{1/2} * C * D^{1/2} * dy, optionally blended with a diagonal
  // (uncorrelated) component via diagonal loading α (default 0.0, i.e. off):
  //   R * dy = (1-α) · D^{1/2} C D^{1/2} · dy + α · D · dy
  // α>0 bounds R's smallest eigenvalue away from zero (faster GMRESR R^-1);
  // diag(R) = D is preserved for any α.
  // where D^{1/2} - diagonal matrix with stddev_ on the diagonal
  //       C - correlations

  // C is estimated using diffusion operator
  // C = Γ * V * W^{-1} * V^{T} * Γ
  // where Γ - normalization operator
  //       V - diffusion operator
  //       W - diagonal matrix with sqrt of row sums of V
  //       V^{T} - adjoint of diffusion operator
  //       W is intrinsically estimated as part of the diffusion operator
  // This is similar to equation (4) in Weaver et al. 2020, QJRMS

  oops::Log::trace() << "ObsErrorDiffusion Multiply: start " << std::endl;

  ASSERT(meshAndNormBuilt_);  // update() must have been called first

  // input copy for the diagonal-loading blend (STEP 6)
  const double alpha = params_.diagonalLoading.value();
  ioda::ObsVector dy_orig(dy);

  // STEP 1: D^{1/2} * dy
  dy *= stddev_;

  // STEP 2: Build a Fieldset and copy dy into it
  // NOTE: This will only work for reading a single variable
  //       (with single channel) from the obsSpace

  const atlas::idx_t nlocs = dy.nlocs();
  const atlas::idx_t npts = geom_->functionSpace().size();

  // localObs2node_ is indexed by local ObsVector location; it must match dy's
  // length (built in update() from the same ObsSpace). Matches randomize()'s guard.
  ASSERT(static_cast<std::size_t>(nlocs) == localObs2node_.size());

  atlas::FieldSet fset;
  for (std::string var : params_.var.value()) {
    atlas::Field obs = geom_->functionSpace().createField<double>(
                       atlas::option::levels(1) | atlas::option::name(var));

    auto obsView = atlas::array::make_view<double, 2>(obs);
    obsView.assign(0.0);
    // Scatter each obs into its mesh node; -1 entries are skipped (control
    // points and ghosts stay 0 and are filled by the diffusion halo exchange).
    for (atlas::idx_t i = 0; i < nlocs; ++i) {
      if (localObs2node_[i] >= 0) obsView(localObs2node_[i], 0) = dy[i];
    }
    // add obs "var" to fieldSet
    fset.add(obs);
  }
  // ============================================================
  // STEP 3: Apply normalization-diffusion-normalization
  //       3.1: Apply normalization Γ (computed in constructor)
  //       3.2: Apply the diffusion operator V and its adjoint V^{T}
  //       3.3: Apply the normalization Γ again
  // Note: W is estimated as part of the diffusion operator
  // ============================================================

  // define helper lambda to apply normalization sqrt
  auto applyNormSqrt = [&](atlas::Field & f) {
    f.haloExchange();
    auto v_f = atlas::array::make_view<double, 2>(f);
    auto v_norm = atlas::array::make_view<double, 2>(this->hzNorm_);
    for (atlas::idx_t i = 0; i < npts; i++) {
      v_f(i, 0) *= v_norm(i, 0);
    }
  };

  // do diffusion i.e. estimate C = Γ * V * W^{-1} * V^{T} * Γ for each field in subset
  for (auto & field : fset) {
    // do diffusion
    applyNormSqrt(field);
    diffusion_->multiplySqrtAD(fset,  oops::Diffusion::Mode::HorizontalOnly);
    diffusion_->multiplySqrtTL(fset,  oops::Diffusion::Mode::HorizontalOnly);
    applyNormSqrt(field);
  }
  // STEP 4: copy data back into dy
  // move this step into diffusion loop (above) for multiple variables?
  for (auto field : fset) {
    auto fieldView = atlas::array::make_view<double, 2>(field);
    // Gather each obs back from its mesh node; -1 entries are inert.
    for (atlas::idx_t i = 0; i < nlocs; ++i) {
      if (localObs2node_[i] >= 0) dy[i] = fieldView(localObs2node_[i], 0);
    }
  }

  // STEP 5: D^{1/2} * C * D^{1/2} * dy
  dy *= stddev_;

  // STEP 6: blend in the diagonal component: dy ← (1-α)·R_corr·dy + α·D·dy
  if (alpha > 0.0) {
    ioda::ObsVector diag_term = dy_orig;
    diag_term *= stddev_;
    diag_term *= stddev_;  // ·stddev_ twice => multiply by variance D
    dy *= (1.0 - alpha);
    dy.axpy(alpha, diag_term);
  }

  oops::Log::trace() << "ObsErrorDiffusion Multiply: end " << std::endl;
}

// -----------------------------------------------------------------------------

void ObsErrorDiffusion::multiply(const ioda::ObsVector & y_in, ioda::ObsVector & y_out) const {
  oops::Log::trace() << "ObsErrorDiffusion multiply(y, yy): starting " << std::endl;

  y_out = y_in;
  multiply(y_out);

  oops::Log::trace() << "ObsErrorDiffusion multiply(y, yy): end " << std::endl;
}

// -----------------------------------------------------------------------------

void ObsErrorDiffusion::inverseMultiply(ioda::ObsVector & dy) const {
  oops::Log::trace() << "ObsErrorDiffusion Inverse Multiply: start " << std::endl;

  ioda::ObsVector dyo = dy;
  dyo.zero();

  int fullInverseIterations_ = params_.InverseIterations.value();
  double fullInverseAccuracy_ = params_.InverseAccuracy.value();

  oops::IdentityMatrix<ioda::ObsVector> Id;

  oops::GMRESR(dyo, dy, *this, Id, fullInverseIterations_, fullInverseAccuracy_);

  dy = dyo;

  oops::Log::trace() << "ObsErrorDiffusion Inverse Multiply: done " << std::endl;
}

// -----------------------------------------------------------------------------

void ObsErrorDiffusion::randomize(ioda::ObsVector & dy) const {
  oops::Log::trace() << "ObsErrorDiffusion Randomize: start " << std::endl;

  ASSERT(meshAndNormBuilt_);  // update() must have been called first

  std::vector<std::string> vars = params_.var.value();

  atlas::FieldSet rand = util::createRandomFieldSet(
        geom_->comm(), geom_->functionSpace(),
        std::vector<size_t>{1}, vars);

  ASSERT(static_cast<std::size_t>(dy.nlocs()) == localObs2node_.size());

  diffusion_->multiplySqrtTL(rand);

  // NOTE:: double check this works for list of more than one var in yaml
  for (std::string var : vars) {
    auto randView = atlas::array::make_view<double, 2>(rand[var]);
    // Gather each obs from its mesh node; -1 entries left untouched.
    for (atlas::idx_t i = 0; i < dy.nlocs(); ++i) {
      if (localObs2node_[i] >= 0) dy[i] = randView(localObs2node_[i], 0);
    }
  }
}

// -----------------------------------------------------------------------------

void ObsErrorDiffusion::save(const std::string & name) const {
  // should this save the normalization coefficients to a file?
  stddev_.save(name);
  oops::Log::trace() << "ObsErrorDiffusion Save: does nothing for now " << std::endl;
}

// -----------------------------------------------------------------------------

std::unique_ptr<ioda::ObsVector> ObsErrorDiffusion::getObsErrors() const {
  return std::make_unique<ioda::ObsVector>(stddev_);
}

// -----------------------------------------------------------------------------

std::unique_ptr<ioda::ObsVector> ObsErrorDiffusion::getInverseVariance() const {
  return std::make_unique<ioda::ObsVector>(inverseVariance_);
}

// -----------------------------------------------------------------------------

void ObsErrorDiffusion::localize(ioda::ObsVector & locvector) const {
  throw eckit::BadParameter("Localizing a correlated R matrix is "
                            "not yet implemented.");
}

// -----------------------------------------------------------------------------

Eigen::MatrixXf ObsErrorDiffusion::localInverseMultiply(const Eigen::MatrixXf & zz) const {
  throw eckit::BadParameter("Localizing a correlated R matrix is "
                            "not yet implemented.");
}

// -----------------------------------------------------------------------------

int ObsErrorDiffusion::localDim() const {
  throw eckit::BadParameter("Localizing a correlated R matrix is "
                            "not yet implemented.");
}

Eigen::VectorXd ObsErrorDiffusion::local_invVarR() const {
  throw eckit::BadParameter("Localizing a correlated R matrix is "
                            "not yet implemented.");
}

// -----------------------------------------------------------------------------
void ObsErrorDiffusion::print(std::ostream & os) const {
  os << "UFO correlated (diffusion) observation error covariance" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsErrorDiffusion::createControlGrid(
    const std::vector<atlas::PointLonLat>& obsnodes,
    const std::vector<int>& obsOwner,
    const double removeWithin,
    const int gridSpacing,
    std::vector<atlas::PointLonLat>& ctrlPoints,
    std::vector<int>& ctrlOwner) const
  {
  oops::Log::trace() << "ObsErrorDiffusion CreateControlGrid: start " << std::endl;
  const atlas::idx_t nObs = obsnodes.size();
  ctrlPoints.clear();
  ctrlOwner.clear();
  //  --- 1. Build control grid ---
  //  --------------------------------
  // number of grid points in longitude
  const int nx = static_cast<int>(std::round(360.0 / gridSpacing));
  // number of grid points in latitude with equator in center
  const int ny = static_cast<int>(std::round(180.0 / gridSpacing)) + 1;
  const atlas::util::Config cfg = atlas::util::Config("type", "regular_lonlat")
        ("nx", nx)
        ("ny", ny)
        ("lon0", 0.0)
        ("lat0", -90.0);
  atlas::RegularLonLatGrid cntrlGrid(cfg);
  //  --- 2. Build KD-tree from obs ---
  //  -------------------------------------
  // KD-tree is built from the GLOBAL obs set, so the control grid is identical on
  // every rank. The payload is the obs index, used both to filter (remove control
  // points near obs) and to assign each kept control point the owning rank of its
  // nearest obs -- keeping the control-point partition geographically consistent.
  atlas::util::IndexKDTree obstree;  //  declare kd tree
  std::vector<size_t> indices(nObs);  //  for payload
  std::iota(indices.begin(), indices.end(), 0);  // assign index to obs payload
  obstree.build(obsnodes, indices);  // build KDTree from obs locations
  //  --- 3. Filter control grid points + inherit owner from nearest obs ---
  //  ----------------------------------------------------------------------
  ctrlPoints.reserve(nx * ny);
  ctrlOwner.reserve(nx * ny);
  int nRemoved = 0;
  for (atlas::idx_t j = 0; j < cntrlGrid.ny(); j++) {
      for (atlas::idx_t i = 0; i < cntrlGrid.nx(); i++) {
          const atlas::PointLonLat gridPt = cntrlGrid.lonlat(i, j);
          const auto result = obstree.closestPoint(gridPt);
          if (result.distance() <= removeWithin) {
              nRemoved++;
          } else {
              ctrlPoints.push_back(gridPt);
              ctrlOwner.push_back(obsOwner[result.payload()]);
          }
      }
  }
  oops::Log::info() << "ObsErrorDiffusion: CreateControlGrid:"
                    << "Control grid total="
                    << nx*ny
                    << " removed=" << nRemoved
                    << " remaining=" << ctrlPoints.size() << std::endl;
  oops::Log::trace() << "ObsErrorDiffusion: CreateControlGrid: end" << std::endl;
  }

}  // namespace ufo
