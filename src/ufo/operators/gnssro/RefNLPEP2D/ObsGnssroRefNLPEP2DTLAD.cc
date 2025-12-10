/*
 * (C) Copyright 2025 Space Sciences and Engineering, LLC (dba PlanetiQ).
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * author Steve Marshall (smarshall@planetiq.com)
 */
#include <cmath>
#include <cstdint>
#include <iterator>
#include <vector>
#include "eckit/exception/Exceptions.h"
#include "ufo/operators/gnssro/RefNLPEP2D/ObsGnssroRefNLPEP2DTLAD.h"
#include "ufo/operators/gnssro/utils/GnssroGeoVaLs.h"
#include "ufo/operators/gnssro/utils/GnssroProfileExtractor.h"
#include "ufo/operators/gnssro/utils/GnssroRayPathGenerator.h"
#include "ufo/utils/Constants.h"
#include "ufo/utils/VertInterp.interface.h"
#include "ufo/variabletransforms/Formulas.h"


#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsGnssroRefNLPEP2DTLAD> makerGnssroRefNLPEP2D_TL_(
    "GnssroRefNLPEP2D");

// -----------------------------------------------------------------------------

ObsGnssroRefNLPEP2DTLAD::ObsGnssroRefNLPEP2DTLAD(const ioda::ObsSpace & odb,
                               const Parameters_ & parameters)
  : LinearObsOperatorBase(odb)
  , varin_()
  , opts_(parameters.options.value())
  , rpParams_(opts_.rayPathGenType, opts_.approxRayLength,
              opts_.top2D, opts_.res, opts_.nHoriz)
  , orchestrator_(obsspace(), rpParams_)
  , refrCalc_(RefractivityCalculator::create(opts_.refrAlgorithm, opts_.useCompress))
  , trajectorySet_(false)
{
  const std::vector<std::string> vv{GnssroGeoVaLs::VAR_TEMP, GnssroGeoVaLs::VAR_SPHUM,
          GnssroGeoVaLs::VAR_PRES, GnssroGeoVaLs::VAR_GPH};
  varin_.reset(new oops::Variables(vv));

  oops::Log::trace() << "ObsGnssroRefNLPEP2DTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsGnssroRefNLPEP2DTLAD::~ObsGnssroRefNLPEP2DTLAD() {
  oops::Log::trace() << "ObsGnssroRefNLPEP2DTLAD destroyed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGnssroRefNLPEP2DTLAD::setTrajectory(
    const GeoVaLs & gv, ObsDiagnostics &, const QCFlags_t &) {
  const char funcName[] = "ObsGnssroRefNLPEP2DTLAD::setTrajectory";
  oops::Log::trace() << funcName << ": started" << std::endl;

  // Return early if there are no geoVals to operate on.
  if (gv.nlocs() == 0)
  {
    oops::Log::trace() << funcName << ": trajectory not set, numGeoVaLs="
                       << gv.nlocs() << std::endl;
    return;
  }

  // Extract the model data from the GeoVaLs.
  GnssroGeoVaLs ggv(gv);

  // Get geometric height and location from the observations.
  std::size_t nlocs = obsspace().nlocs();
  std::vector<float> latitudes(nlocs);
  std::vector<double> geomHeights(nlocs);
  obsspace().get_db("MetaData", "latitude", latitudes);
  obsspace().get_db("MetaData", "height", geomHeights);

  std::size_t profileIdx = 0;
  std::size_t ovecIdx = 0;  // Index into locations in ObsSpace (nlocs)
  double missing = util::missingValue<double>();

  //  Iterate over RO profiles.
  const GnssroRayPathOrchestrator::Profiles_& profiles = orchestrator_.profiles();
  for (const auto& roProfile : profiles)
  {
    int seqNum = roProfile->seqNum();
    oops::Log::debug() << funcName << ": processing profile seqNum=" << seqNum << std::endl;

    // Iterate over rays within this profile.
    for (std::size_t rayIdx = 0; rayIdx < roProfile->numRays(); ++rayIdx, ++ovecIdx)
    {
      const GnssroRayPathOrchestrator::Ray_& ray = roProfile->ray(rayIdx);

      std::size_t tpNodeIdx = ray->tpNodeIdx();

      // Initialize the trajectories for this ray.
      ray->traj().initialize(ray->numNodes());

      for (std::size_t nodeIdx = 0; nodeIdx < ray->numNodes(); ++nodeIdx)
      {
        // Include segment length and unit conversion in factor for refractivity Jacobian.
        double segmentLength = ray->segLen(nodeIdx) * RefractivityCalculator::N_TO_EXCESS_IOR;

        // Get geometric height of this node
        double geomHgt;
        if (nodeIdx == tpNodeIdx) {
            geomHgt = geomHeights[ovecIdx];
        } else {  // Average of edges of segment.
            double thisGeomHgt = ray->geomHgt(nodeIdx);
            double nextGeomHgt = ray->geomHgt(nodeIdx + 1);
            geomHgt = 0.5 * (thisGeomHgt + nextGeomHgt);
        }

        // Get trajectory information from the model profile.
        GnssroRayPathOrchestrator::TrajTuple trajTuple =
            GnssroRayPathOrchestrator::modelTrajectory(
                geomHgt, latitudes[ovecIdx], ggv, profileIdx, refrCalc_);

        // Save the weighted Jacobian values.
        ray->traj().jacobianT(nodeIdx) = segmentLength * trajTuple.dNdT();
        ray->traj().jacobianP(nodeIdx) = segmentLength * trajTuple.dNdP();
        ray->traj().jacobianQ(nodeIdx) = segmentLength * trajTuple.dNdQ();

        // Save the vertical interpolation weights in the trajectory.
        ray->traj().wf(nodeIdx) = trajTuple.wf();
        ray->traj().wi(nodeIdx) = trajTuple.wi();
        oops::Log::debug() << funcName << ": seqNum=" << seqNum << ", rayIdx=" << rayIdx
                  << ", nodeIdx=" << nodeIdx << " (" << tpNodeIdx << "): 2D Traj"
                  << ", jacT=" << ray->traj().jacobianT(nodeIdx)
                  << ", jacQ=" << ray->traj().jacobianQ(nodeIdx)
                  << ", jacP=" << ray->traj().jacobianP(nodeIdx)
                  << ", wf=" << ray->traj().wf(nodeIdx)
                  << ", wi=" << ray->traj().wi(nodeIdx) << std::endl;

        // If using spherical symmetry, all nodes use the same model column.
        if (!ray->assumeSphericalSymmetry()) {
          ++profileIdx;
        }
      }  // end of iteration of nodes in a single ray.

      // Mark trajectory for this ray as complete and go to next observed location.
      ray->traj().finalize();
    }  // end of iteration of rays in a single RO profile

    oops::Log::debug() << funcName << ": seqNum=" << seqNum << ", preparing for next profile"
          << ", ovecIdx=" << ovecIdx << std::endl;
  }  // end of iteration over RO profiles

  trajectorySet_ = true;
  oops::Log::trace() << funcName << ": trajectory set" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGnssroRefNLPEP2DTLAD::simulateObsTL(
    const GeoVaLs & gv, ioda::ObsVector & ovec, const QCFlags_t & qc_flags) const
{
  const char funcName[] = "ObsGnssroRefNLPEP2DTLAD::simulateObsTL";
  oops::Log::trace() << funcName << ": started" << std::endl;

  // Return early if there are no geoVals to operate on.
  if (gv.nlocs() == 0)
  {
    oops::Log::trace() << funcName << ": not run because numGeoVaLs="
                       << gv.nlocs() << std::endl;
    return;
  }

  //  Validate our state.
  if (!trajectorySet_)
  {
    oops::Log::error() << funcName << ": trajectory was not set, nlocs=" << gv.nlocs()
                       << std::endl;
    throw eckit::SeriousBug("ObsGnssroRefNLPEP2DTLAD::simulateObsTL: trajectory was not set",
                            Here());
  }

  if (gv.nlocs() != ovec.size())
  {
    oops::Log::error() << funcName << ": GeoVaLs:nlocs=" << gv.nlocs() << " not equal to "
                       << "size of ObsVector=" << ovec.size() << std::endl;
    throw eckit::SeriousBug("ObsGnssroRefNLPEP2DTLAD::simulateObsTL: nlocs in geovals "
                            "not equal to size of ObsVector", Here());
  }

  // Extract the model data from the GeoVaLs.
  GnssroGeoVaLs ggv(gv);

  const double missing = util::missingValue<double>();
  std::size_t profileIdx = 0;
  std::size_t ovecIdx = 0;  // Index into locations in ObsVector (nlocs)

  //  Iterate over RO profiles.
  const GnssroRayPathOrchestrator::Profiles_& profiles = orchestrator_.profiles();
  for (const auto& roProfile : profiles)
  {
    int seqNum = roProfile->seqNum();
    oops::Log::debug() << funcName << ": processing profile seqNum=" << seqNum << std::endl;

    // Determine range of observation indices in this profile.
    std::size_t startOvecIdx = ovecIdx;
    std::size_t endOvecIdx = ovecIdx + roProfile->numRays();
    ovecIdx = endOvecIdx;

    // Get range of model profile indices in this RO profile
    std::size_t startProfileIdx = profileIdx;
    std::size_t endProfileIdx = startProfileIdx + roProfile->totalSampledLocations();
    profileIdx = endProfileIdx;

    oops::Log::debug() << funcName << ": seqNum=" << seqNum << ": processing profile "
            << "with numRays=" << roProfile->numRays() << " ovecIdx range=["
            << startOvecIdx << ":" << endOvecIdx << "), profile-wide profileIdx range=["
            << startProfileIdx << "," << endProfileIdx << "), starting at ovecIdx="
            << ovecIdx - 1 << " and proceeding backwards" << std::endl;

    // Iterate over rays within this profile from top to bottom (reverse order)
    for (GnssroProfileRayPaths::const_reverse_iterator rayItr = roProfile->crbegin();
         rayItr != roProfile->crend(); ++rayItr)
    {
      // Get the integer rayIdx for use in diagnostic messages.
      // Get a forward iterator pointing to one position after the reverse iterator.
      // Then convert the forward iterator into an integer position.
      // This prevents us from even computing a negative rayIdx, which would happen
      // if we used an integer index as a reverse iterator.
      --ovecIdx;
      auto fwdItr = rayItr.base() - 1;
      std::size_t rayIdx = std::distance(roProfile->cbegin(), fwdItr);
      const GnssroProfileRayPaths::Ray_& ray = *rayItr;

      // Get the model profile range for this ray
      std::size_t rayEndProfileIdx = profileIdx;
      std::size_t rayStartProfileIdx = rayEndProfileIdx - ray->numLocations();
      profileIdx = rayStartProfileIdx;
      oops::Log::debug() << funcName << ": seqNum=" << seqNum << ", rayIdx="
            << rayIdx << ", numNodes=" << ray->numNodes() << ", numSampLocs="
            << ray->numLocations() << ", ray profileIdx range=[" << rayStartProfileIdx
            << "," << rayEndProfileIdx << "), starting at profileIdx="
            << profileIdx << std::endl;

      std::size_t tpNodeIdx = ray->tpNodeIdx();

      double fwdTL = 0.0;
      ovec[ovecIdx] = 0.0;  // Accumulate contribution to TL hofx from each node.
      for (std::size_t nodeIdx = 0; nodeIdx < ray->numNodes(); ++nodeIdx)
      {
        // Get references to model vertical profiles (columns) for this profileIdx.
        const GnssroGeoVaLs::Profile_& tempProf = ggv.tempProfile(profileIdx);
        const GnssroGeoVaLs::Profile_& sphumProf = ggv.sphumProfile(profileIdx);
        const GnssroGeoVaLs::Profile_& presProf = ggv.presProfile(profileIdx);

        // Vertically interpolate using previously computed weights for this node.
        int wi = ray->traj().wi(nodeIdx);
        double wf = ray->traj().wf(nodeIdx);
        double temp;    // Interpolated temperature perturbation in K
        double sphum;   // Interpolated specific humidity perturbation in kg/kg
        double pres;    // Interpolated pressure perturbation in Pa
        vert_interp_apply_f90(ggv.nlevs(), tempProf.data(), temp, wi, wf);
        vert_interp_apply_f90(ggv.nlevs(), sphumProf.data(), sphum, wi, wf);
        vert_interp_apply_f90(ggv.nlevs(), presProf.data(), pres, wi, wf);

        // Compute the forward part
        // Note: Jacobian at nodes is already node-weighted.
        double fwdIncr = (
                  ray->traj().jacobianT(nodeIdx) * temp
                + ray->traj().jacobianP(nodeIdx) * pres
                + ray->traj().jacobianQ(nodeIdx) * sphum);
        fwdTL += fwdIncr;

        // If using spherical symmetry, all nodes use the same model column.
        if (!ray->assumeSphericalSymmetry()) {
          ++profileIdx;
        }
      }  // end of loop over nodes

      ovec[ovecIdx] = fwdTL;

      // Prepare for the next lower ray by setting the profileIdx to point to
      // the start of the model profiles for the current ray, which is also the
      // end of the range for the next lower ray within the profile.
      profileIdx = rayStartProfileIdx;
    }  // end of iteration of rays in a single RO profile

    // Prepare for iteration into next RO profile.
    oops::Log::debug() << funcName << ": seqNum=" << seqNum << " DONE: "
          << "ovecIdx=" << ovecIdx << std::endl;
    ovecIdx = endOvecIdx;
    if (profileIdx != startProfileIdx) {
      oops::Log::error() << funcName << ": seqNum=" << seqNum
              << ": model profile indexing error: after RO profile processed,"
              << " profileIdx=" << profileIdx << "; should be same as "
              << "startProfileIdx=" << startProfileIdx << std::endl;
    }
    profileIdx = endProfileIdx;
    oops::Log::debug() << funcName << ": seqNum=" << seqNum << " DONE: "
              << "skipping to next profile ovecIdx=" << ovecIdx
              << ", profileIdx=" << profileIdx << std::endl;
  }  // end of iteration over RO profiles

  oops::Log::trace() << funcName << ": complete" << std::endl;
}

//----------------------------------------------------------------------------

void ObsGnssroRefNLPEP2DTLAD::simulateObsAD(GeoVaLs & gv, const ioda::ObsVector & ovec,
                                   const QCFlags_t & qc_flags) const
{
  const char funcName[] = "ObsGnssroRefNLPEP2DTLAD::simulateObsAD";
  oops::Log::trace() << funcName << ": started" << std::endl;

  // Return early if there are no geoVals to operate on.
  if (gv.nlocs() == 0)
  {
    oops::Log::trace() << funcName << ": not run because numGeoVaLs="
                       << gv.nlocs() << std::endl;
    return;
  }

  //  Validate our state.
  if (!trajectorySet_)
  {
    oops::Log::error() << funcName << ": trajectory was not set, nlocs=" << gv.nlocs()
                       << std::endl;
    throw eckit::SeriousBug("ObsGnssroRefNLPEP2DTLAD::simulateObsTL: trajectory was not set",
                            Here());
  }

  if (gv.nlocs() != ovec.size())
  {
    oops::Log::error() << funcName << ": GeoVaLs:nlocs=" << gv.nlocs() << " not equal to "
                       << "size of ObsVector=" << ovec.size() << std::endl;
    throw eckit::SeriousBug("ObsGnssroRefNLPEP2DTLAD::simulateObsTL: nlocs in geovals "
                            "not equal to size of ObsVector", Here());
  }

  // Extract the model data from the GeoVaLs.
  GnssroGeoVaLs ggv(gv);

  const double missing = util::missingValue<double>();
  std::size_t profileIdx = 0;
  std::size_t ovecIdx = 0;  // Index into locations in ObsVector (nlocs)
  //  Iterate over RO profiles.
  const GnssroRayPathOrchestrator::Profiles_& profiles = orchestrator_.profiles();
  for (const auto& roProfile : profiles)
  {
    int seqNum = roProfile->seqNum();
    oops::Log::debug() << funcName << ": processing profile seqNum=" << seqNum << std::endl;

    // Determine range of observation indices in this profile.
    std::size_t startOvecIdx = ovecIdx;
    std::size_t endOvecIdx = ovecIdx + roProfile->numRays();
    ovecIdx = endOvecIdx;

    // Get range of model profile indices in this RO profile
    std::size_t startProfileIdx = profileIdx;
    std::size_t endProfileIdx = startProfileIdx + roProfile->totalSampledLocations();
    profileIdx = endProfileIdx;

    oops::Log::debug() << funcName << ": seqNum=" << seqNum << ": processing profile "
            << "with numRays=" << roProfile->numRays() << " ovecIdx range=["
            << startOvecIdx << ":" << endOvecIdx << "), profile-wide profileIdx range=["
            << startProfileIdx << "," << endProfileIdx << "), starting at ovecIdx="
            << ovecIdx - 1 << " and proceeding backwards" << std::endl;

    // Iterate over rays within this profile from top to bottom (reverse order)
    for (GnssroProfileRayPaths::const_reverse_iterator rayItr = roProfile->crbegin();
         rayItr != roProfile->crend(); ++rayItr)
    {
      // Get the integer rayIdx for use in diagnostic messages.
      // Get a forward iterator pointing to one position after the reverse iterator.
      // Then convert the forward iterator into an integer position.
      // This prevents us from even computing a negative rayIdx, which would happen
      // if we used an integer index as a reverse iterator.
      --ovecIdx;
      auto fwdItr = rayItr.base() - 1;  // fwdItr points to same loc as reverse itr
      std::size_t rayIdx = std::distance(roProfile->cbegin(), fwdItr);
      const GnssroProfileRayPaths::Ray_& ray = *rayItr;

      // Get the model profile range for this ray
      std::size_t rayEndProfileIdx = profileIdx;
      std::size_t rayStartProfileIdx = rayEndProfileIdx - ray->numLocations();
      profileIdx = rayStartProfileIdx;
      oops::Log::debug() << funcName << ": seqNum=" << seqNum << ", rayIdx="
            << rayIdx << ", ovecIdx=" << ovecIdx << ", ovec.size=" << ovec.size()
              << ", numNodes=" << ray->numNodes() << ", numSampLocs="
              << ray->numLocations() << ", ray profileIdx range=[" << rayStartProfileIdx
              << "," << rayEndProfileIdx << "), starting at profileIdx="
              << profileIdx << std::endl;

      std::size_t tpNodeIdx = ray->tpNodeIdx();

      for (std::size_t nodeIdx = 0; nodeIdx < ray->numNodes(); ++nodeIdx)
      {
        // Compute the perturbations in each variable for the forward part
        double dTemp = ovec[ovecIdx] * ray->traj().jacobianT(nodeIdx);
        double dSphum = ovec[ovecIdx] * ray->traj().jacobianQ(nodeIdx);
        double dPres = ovec[ovecIdx] * ray->traj().jacobianP(nodeIdx);

        // These perturbations should be divided over one or two vertical
        // levels in the profile (wi and wi+1), with relative amounts
        // determined by the wf value. This logic follows the Fortran
        // function vert_interp_apply_ad.
        // NOTE: when weighting the contributions, we want to include both
        // the horizontal weight (segmentLength) and the vertical weight (wf).
        // However, segmentLength is already part of the value returned
        // by the node-specific Jacobians.
        int wi = ray->traj().wi(nodeIdx);
        int cwi = wi - 1;  // 0-based index equivalent to wi.
        double wf = ray->traj().wf(nodeIdx);
        if (cwi >= 0)  // Do not go below lowest vertical level in model
        {
          std::size_t levelIdx = cwi;
          double weight = wf;
          oops::Log::debug() << funcName << ": seqNum=" << seqNum
                    << ", rayIdx=" << rayIdx << ", nodeIdx=" << nodeIdx
                    << ", prof=" << profileIdx << ", lev=" << levelIdx
                    << ", incr temp=" << ggv.temp(profileIdx, levelIdx)
                    << " by " << (weight * dTemp) << ", incr sphum="
                    << ggv.sphum(profileIdx, levelIdx) << " by "
                    << (weight * dSphum) << ", incr pres="
                    << ggv.pres(profileIdx, levelIdx) << " by "
                    << (weight * dPres) << std::endl;
          double & geoval_temp = ggv.temp(profileIdx, levelIdx);
          if (geoval_temp != missing) {
            geoval_temp += (weight * dTemp);
          }
          double & geoval_sphum = ggv.sphum(profileIdx, levelIdx);
          if (geoval_sphum != missing) {
            geoval_sphum += (weight * dSphum);
          }
          double & geoval_pres = ggv.pres(profileIdx, levelIdx);
          if (geoval_pres != missing) {
            geoval_pres += (weight * dPres);
          }
          oops::Log::debug() << funcName << ": seqNum=" << seqNum
                    << ", rayIdx=" << rayIdx << ", nodeIdx=" << nodeIdx
                    << ", prof=" << profileIdx << ", lev=" << levelIdx
                    << ", set temp=" << ggv.temp(profileIdx, levelIdx)
                    << ", set sphum=" << ggv.sphum(profileIdx, levelIdx)
                    << ", set pres=" << ggv.pres(profileIdx, levelIdx)
                    << std::endl;
        }

        if (cwi + 1 < ggv.nlevs())  // Do not go above highest level in model
        {
          std::size_t levelIdx = cwi + 1;
          double weight = (1.0 - wf);
          oops::Log::debug() << funcName << ": seqNum=" << seqNum
                    << ", rayIdx=" << rayIdx << ", nodeIdx=" << nodeIdx
                    << ", prof=" << profileIdx << ", lev=" << levelIdx
                    << ", incr temp=" << ggv.temp(profileIdx, levelIdx)
                    << " by " << (weight * dTemp) << ", incr sphum="
                    << ggv.sphum(profileIdx, levelIdx) << " by "
                    << (weight * dSphum) << ", incr pres="
                    << ggv.pres(profileIdx, levelIdx) << " by "
                    << (weight * dPres) << std::endl;
          double & geoval_temp = ggv.temp(profileIdx, levelIdx);
          if (geoval_temp != missing) {
            geoval_temp += (weight * dTemp);
          }
          double & geoval_sphum = ggv.sphum(profileIdx, levelIdx);
          if (geoval_sphum != missing) {
            geoval_sphum += (weight * dSphum);
          }
          double & geoval_pres = ggv.pres(profileIdx, levelIdx);
          if (geoval_pres != missing) {
            geoval_pres += (weight * dPres);
          }
          oops::Log::debug() << funcName << ": seqNum=" << seqNum
                    << ", rayIdx=" << rayIdx << ", nodeIdx=" << nodeIdx
                    << ", prof=" << profileIdx << ", lev=" << levelIdx
                    << ", set temp=" << ggv.temp(profileIdx, levelIdx)
                    << ", set sphum=" << ggv.sphum(profileIdx, levelIdx)
                    << ", set pres=" << ggv.pres(profileIdx, levelIdx)
                    << std::endl;
        }

        // If using spherical symmetry, all nodes use the same model column.
        if (!ray->assumeSphericalSymmetry()) {
          ++profileIdx;
        }
      }  // end of iteration over nodes for adjoint.

      // Prepare for the next lower ray by setting the profileIdx to point to
      // the start of the model profiles for the current ray, which is also the
      // end of the range for the next lower ray within the profile.
      profileIdx = rayStartProfileIdx;
    }  // end of iteration of rays in a single RO profile

    //  Prepare for iteration into next RO profile.
    oops::Log::debug() << funcName << ": seqNum=" << seqNum << " DONE: "
          << "ovecIdx=" << ovecIdx << std::endl;
    ovecIdx = endOvecIdx;
    if (profileIdx != startProfileIdx) {
      oops::Log::error() << funcName << ": seqNum=" << seqNum
              << ": model profile indexing error: after RO profile processed,"
              << " profileIdx=" << profileIdx << "; should be same as "
              << "startProfileIdx=" << startProfileIdx << std::endl;
    }
    profileIdx = endProfileIdx;
    oops::Log::debug() << funcName << ": seqNum=" << seqNum << " DONE: "
              << "skipping to next profile ovecIdx=" << ovecIdx
              << ", profileIdx=" << profileIdx << std::endl;
  }  // end of iteration over RO profiles

  // Save the data from the GnssroGeoVaLs to the GeoVaLs passed to us.
  ggv.saveTo(gv);

  oops::Log::trace() << funcName << ": complete" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGnssroRefNLPEP2DTLAD::print(std::ostream & os) const {
  os << "ObsGnssroRefNLPEP2DTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
