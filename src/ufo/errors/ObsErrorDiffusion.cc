
/*
 * (C) Copyright 2025 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <utility>
#include <vector>

#include "atlas/grid/Grid.h"
#include "atlas/mesh/Mesh.h"
#include "atlas/mesh/MeshBuilder.h"
#include "atlas/output/Gmsh.h"
#include "atlas/util/Geometry.h"
#include "atlas/util/KDTree.h"

#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"

#include "oops/assimilation/GMRESR.h"
#include "oops/base/FieldSet3D.h"
#include "oops/base/IdentityMatrix.h"
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
  ASSERT(comm_.size() == 1);  // does not support mesh generation on multiple PEs
  inverseVariance_ = stddev_;
  inverseVariance_ *= stddev_;
  inverseVariance_.invert();
  oops::Log::trace() << "ObsErrorDiffusion::ObsErrorDiffusion constructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsErrorDiffusion::update(const ioda::ObsVector & obserr) {
  oops::Log::trace() << "ObsErrorDiffusion update: start " << std::endl;

  stddev_ = obserr;
  inverseVariance_ = stddev_;
  inverseVariance_ *= stddev_;
  inverseVariance_.invert();

  // *for single channel/var* all obs locations on this PE (including obs not passing qc)
  const int nlocs = obserr.nlocs();
  // *for single channel/var* nobs is number of obs passing QC across all mpi tasks
  const int nobs = obserr.nobs();
  std::vector<float> lons(nlocs);
  std::vector<float> lats(nlocs);
  std::vector<atlas::PointLonLat> obsnodes;
  obsnodes.reserve(nlocs);
  std::vector<double> gridXY;  // this will have length of 2*(nobs + control_grid_size)
  gridXY.reserve(2*nlocs);  // overestimate with QC; underestimate (possibly) with a control grid

  // lons/lats arrays contains ALL obs locations
  obserr.space().get_db("MetaData", "longitude", lons);
  obserr.space().get_db("MetaData", "latitude", lats);

  const double missing = util::missingValue<double>();

  // loop over all the obs
  for (int i = 0; i < nlocs; ++i) {
    if (obserr[i] != missing) {  // currently, passivated obs are set to missingValue
      // fill vector in form: [ lon_0, lat_0, lon_1, lat_1, ... ]
      gridXY.push_back(lons[i]);
      gridXY.push_back(lats[i]);
      atlas::PointLonLat p(lons[i], lats[i]);
      p.normalise();
      obsnodes.push_back(p);
    }
  }

  // build config to create mesh
  eckit::LocalConfiguration fspaceConfig;
  fspaceConfig.set("function space", "NodeColumns");  // always will be NodeColumns
  fspaceConfig.set("grid.type", "unstructured");      // always will be unstructured
  fspaceConfig.set("grid.xy", gridXY);
  fspaceConfig.set("partitioner", "equal_regions");
  fspaceConfig.set("no point on last task", true);

  atlas::Grid grid{};
  atlas::FunctionSpace fspace{};
  atlas::grid::Partitioner partitioner{};
  atlas::Mesh mesh{};
  atlas::FieldSet fset{};

  oops::Log::trace() << "Constructing diffusion mesh with " << obserr.nobs()
                     << " valid observations" << std::endl;

  //------------------------------------------------------------------------------
  // If control grid parameters are provided,
  // use them to create a coarser grid for the diffusion operator
  //
  // NOTE: not setting control grid as default as the grid spacing and
  // remove within parameters should be explicitly set,
  // depending on the observation network geometry,
  // to ensure the control grid is created as intended
  //------------------------------------------------------------------------------
  if (params_.controlGrid.value()) {
    const int gridSpacing = params_.controlGrid.value()->gridSpacing.value();
    const double removeWithin = params_.controlGrid.value()->removeWithin.value();

    // returns merged obs + remaining control grid points
    std::vector<atlas::PointLonLat> mergedPoints =
               createControlGrid(obsnodes, removeWithin, gridSpacing);
    // nMerged is (final number of control points) + (nobs)
    const atlas::idx_t nMerged = mergedPoints.size();
    // obsOffset_ is number of control points
    obsOffset_ = nMerged - nobs;

    // update gridXY with merged points for mesh creation
    gridXY.resize(2 * mergedPoints.size());
    for (atlas::idx_t k = 0; k < nMerged; k++) {
        gridXY[2*k]   = mergedPoints[k].lon();
        gridXY[2*k+1] = mergedPoints[k].lat();
    }
    fspaceConfig.set("grid.xy", gridXY);
  }
  // Uses Delauney triangulation for mesh generation
  util::setupFunctionSpace(comm_, fspaceConfig, grid, partitioner, mesh, fspace, fset);
  geom_.reset(new oops::GeometryData(fspace, fset, true, comm_));
  diffusionGeom_ = oops::Diffusion::calculateDerivedGeom(std::as_const(*geom_));
  diffusion_.reset(new oops::Diffusion(std::as_const(*geom_), diffusionGeom_));

  const atlas::idx_t npts = geom_->functionSpace().size();

  // Save mesh and point-type marker field for visualization with gmsh
  if (params_.outputDiffusionMesh.value()) {
    // Create marker field: 0 = control grid point, 1 = obs point
    atlas::Field markerField = geom_->functionSpace().createField<double>(
        atlas::option::levels(1) | atlas::option::name("pointType"));
    auto v_marker = atlas::array::make_view<double, 2>(markerField);
    v_marker.assign(0.0);

    for (atlas::idx_t i = 0; i < nobs; ++i) {
        // control grid points at BEGINNING of array
        v_marker(obsOffset_ + i, 0) = 1.0;
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
  const double val = params_.lscale;
  v_hzScales.assign(val);
  // Set the scales for control points to 0 so they don't influence the diffusion of obs points
  for (atlas::idx_t i = 0; i < obsOffset_; ++i) {
    v_hzScales(i, 0) = 0;
  }
  atlas::FieldSet scales;
  scales.add(hzScales);

  oops::Log::info() << "  hzScales created with fixed Gaussian sigma = "
                    << val << std::endl;
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
  oops::Log::trace() << "ObsErrorDiffusion update: end" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsErrorDiffusion::multiply(ioda::ObsVector & dy) const {
  // R * dy = D^{1/2} * C * D^{1/2} * dy
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

  // STEP 1: D^{1/2} * dy
  dy *= stddev_;

  // STEP 2: Build a Fieldset and copy dy into it
  // NOTE: This will only work for reading a single variable
  //       (with single channel) from the obsSpace

  const atlas::idx_t nobs = dy.nobs();  // store as class member? works with mult vars??
  const atlas::idx_t nlocs = dy.nlocs();
  const atlas::idx_t npts = geom_->functionSpace().size();
  const double missing = util::missingValue<double>();

  atlas::FieldSet fset;
  for (std::string var : params_.var.value()) {
    atlas::Field obs = geom_->functionSpace().createField<double>(
                       atlas::option::levels(1) | atlas::option::name(var));

    auto obsView = atlas::array::make_view<double, 2>(obs);
    obsView.assign(0.0);
    atlas::idx_t validObsIn = 0;
    // dy has length = nlocs (for single var/channel)
    for (atlas::idx_t i = 0; i < nlocs; ++i) {
      if (dy[i] != missing) {
        // only copy obs from valid locations
        obsView(obsOffset_ + validObsIn, 0) = dy.data()[i];
        ++validObsIn;
      }
    }
    ASSERT(validObsIn == nobs);
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
    atlas::idx_t validObsOut = 0;
    // using nlocs to fully iterate through dy
    for (atlas::idx_t i = 0; i < nlocs; ++i) {
      if (dy[i] != missing) {
        dy[i] = fieldView(obsOffset_ + validObsOut, 0);
        ++validObsOut;
      }
    }
    ASSERT(validObsOut == nobs);
  }

  // STEP 5: D^{1/2} * C * D^{1/2} * dy
  dy *= stddev_;

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

  int fullInverseIterations_ = params_.InverseIterations.value();  //, 10);
  double fullInverseAccuracy_ = params_.InverseAccuracy.value();   //, 1.0e-5);

  oops::IdentityMatrix<ioda::ObsVector> Id;

  oops::GMRESR(dyo, dy, *this, Id, fullInverseIterations_, fullInverseAccuracy_);

  dy = dyo;

  oops::Log::trace() << "ObsErrorDiffusion Inverse Multiply: done " << std::endl;
}

// -----------------------------------------------------------------------------

void ObsErrorDiffusion::randomize(ioda::ObsVector & dy) const {
  oops::Log::trace() << "ObsErrorDiffusion Randomize: start " << std::endl;

  std::vector<std::string> vars = params_.var.value();

  atlas::FieldSet rand = util::createRandomFieldSet(
        geom_->comm(), geom_->functionSpace(),
        std::vector<size_t>{1}, vars);

  const atlas::idx_t npts = geom_->functionSpace().size();  // can do this with rand?
  const atlas::idx_t nobs = npts - obsOffset_;

  ASSERT(nobs == dy.nobs());  // ensure valid obs locs in diffusion grid equals num_valid obs in dy
  ASSERT(nobs == dy.nlocs());  // assumes no QC has been performed (all obs locations are valid)

  diffusion_->multiplySqrtTL(rand);

  // NOTE:: double check this works for list of more than one var in yaml
  for (std::string var : vars) {
    auto randView = atlas::array::make_view<double, 2>(rand[var]);
    // iterating over nlocs will not work if obs have been QC'ed out of dy, but
    // since this is only used in the obsErrorCovariance test, this should be okay
    for (atlas::idx_t i = 0; i < dy.nlocs(); ++i) {
      dy[i] = randView(obsOffset_ + i, 0);
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
  os << "UFO Diagonal observation error covariance, inverse variances: " << std::endl;
}

// -----------------------------------------------------------------------------

std::vector <atlas::PointLonLat> ObsErrorDiffusion::createControlGrid(
    const std::vector<atlas::PointLonLat>& obsnodes,
    const double removeWithin,
    const int gridSpacing)
  {
  oops::Log::trace() << "ObsErrorDiffusion CreateControlGrid: start " << std::endl;
  const atlas::idx_t nlocs = obsnodes.size();
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
  atlas::util::IndexKDTree obstree;  //  declare kd tree
  std::vector<size_t> indices(nlocs);  //  for payload
  std::iota(indices.begin(), indices.end(), 0);  // assign index to obs payload
  obstree.build(obsnodes, indices);  // build KDTree from obs locations
  //  --- 3. Filter control grid points ---
  //  --------------------------------------
  std::vector<atlas::PointLonLat> remainingGridPoints;
  remainingGridPoints.reserve(nx * ny);
  int nRemoved = 0;
  for (atlas::idx_t j = 0; j < cntrlGrid.ny(); j++) {
      for (atlas::idx_t i = 0; i < cntrlGrid.nx(); i++) {
          const atlas::PointLonLat gridPt = cntrlGrid.lonlat(i, j);
          const auto result = obstree.closestPoint(gridPt);
          if (result.distance() <= removeWithin) {
              nRemoved++;
              oops::Log::debug() << "ObsErrorDiffusion: CreateControlGrid:"
                                 << "Removed grid point ("
                                 << gridPt.lon() << ", " << gridPt.lat()
                                 << ") — nearest obs #" << result.payload()
                                 << " (lon=" << obsnodes[result.payload()].lon()
                                 << ", lat=" << obsnodes[result.payload()].lat()
                                 << ") at distance=" << result.distance() << "m\n";
          } else {
              remainingGridPoints.push_back(gridPt);
          }
      }
  }
  oops::Log::info() << "ObsErrorDiffusion: CreateControlGrid:"
                    << "Control grid total="
                    << nx*ny
                    << " removed=" << nRemoved
                    << " remaining=" << remainingGridPoints.size() << std::endl;
  // --- 4. Merge remaining control grid points with obs locations ---
  //---------------------------------------------------------------------
  std::vector<atlas::PointLonLat> mergedPoints;
  mergedPoints.reserve(remainingGridPoints.size() + nlocs);
  mergedPoints.insert(mergedPoints.end(),
                    remainingGridPoints.begin(), remainingGridPoints.end());
  mergedPoints.insert(mergedPoints.end(),
                    obsnodes.begin(), obsnodes.end());
  oops::Log::trace() << "ObsErrorDiffusion: CreateControlGrid: end" << std::endl;
  return mergedPoints;
  }

}  // namespace ufo
