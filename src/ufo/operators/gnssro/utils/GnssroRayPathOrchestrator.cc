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
#include <utility>
#include "ioda/ObsSpace.h"
#include "ufo/operators/gnssro/utils/GnssroGeoVaLs.h"
#include "ufo/operators/gnssro/utils/GnssroRayPathOrchestrator.h"
#include "ufo/operators/gnssro/utils/RefractivityCalculator.h"
#include "ufo/utils/Constants.h"
#include "ufo/utils/VertInterp.interface.h"
#include "ufo/variabletransforms/Formulas.h"

namespace ufo {

// -----------------------------------------------------------------------------
static const char className[] = "GnssroRayPathOrchestrator";

GnssroRayPathOrchestrator::GnssroRayPathOrchestrator(
        const ioda::ObsSpace & odb,
        const GnssroRayPathParameters & params)
  : odb_(odb)
  , params_(params)
  , extractor_(odb_)
  , generator_()
  , profiles_()
{
  oops::Log::trace() << className << " ctor enter" << std::endl;

  oops::Log::debug() << className << ": " << extractor_.slices().size()
                     << " profiles found" << std::endl;

  if (extractor_.is_empty()) {
      return;  // No observations to process.
  }

  // Construct the ray path generator.
  generator_ = GnssroRayPathGenerator::create(odb, params_);

  // Create a set of ray paths for each profile.
  std::size_t numRays = 0;
  for (GnssroProfileExtractor::const_iterator itr = extractor_.cbegin();
       itr != extractor_.cend(); ++itr)
  {
    oops::Log::debug() << className << ": Profile slice "
                       << itr->seqNum() << " from " << itr->start()
                       << " to " << itr->end() << std::endl;
    Profile_ profileRayPaths = generator_->makeProfileRayPaths(*itr);
    numRays += profileRayPaths->numRays();
    profiles_.push_back(std::move(profileRayPaths));
  }

  oops::Log::debug() << className << ": generated " << numRays << " rays for "
                     << profiles_.size() << " profiles" << std::endl;
  oops::Log::trace() << className << " ctor exit" << std::endl;
}
// -----------------------------------------------------------------------------

GnssroRayPathOrchestrator::~GnssroRayPathOrchestrator()
{
  oops::Log::trace() << className << " dtor called" << std::endl;
}
// -----------------------------------------------------------------------------

std::size_t GnssroRayPathOrchestrator::numSampledLocs() const
{
  // We will have one sampled location for each node.
  // Iterate over RO profiles
  std::size_t numLocs = 0;
  for (Profiles_::const_iterator p = profiles_.cbegin(); p != profiles_.cend(); ++p)
  {
      // Iterate over rays in a profile
      const Profile_& roProfile = *p;
      for (std::size_t rayIdx = 0; rayIdx < roProfile->numRays(); ++rayIdx)
      {
          numLocs += roProfile->ray(rayIdx)->numLocations();
      }
  }
  return numLocs;
}
// -----------------------------------------------------------------------------

std::size_t GnssroRayPathOrchestrator::getSampledLatLonTimes(
    std::vector<float>& lats,
    std::vector<float>& lons,
    std::vector<util::DateTime>& times) const
{
  const char funcName[] = "GnssroRayPathOrchestrator::getSampledLatLonTimes";
  lats.clear();
  lons.clear();
  times.clear();
  std::size_t numLocs = numSampledLocs();
  lats.reserve(numLocs);
  lons.reserve(numLocs);
  times.reserve(numLocs);

  oops::Log::debug() << funcName << ": Number of profiles in extractor="
          << extractor_.slices().size() << ", and number of profiles in profiles= "
          << profiles_.size()  << std::endl;

  GnssroProfileExtractor::const_iterator extItr = extractor_.cbegin();
  for (Profiles_::const_iterator itr = profiles_.cbegin(); itr != profiles_.cend(); ++itr)
  {
    const Profile_& roProfile = *itr;
    for (std::size_t rayIdx = 0; rayIdx < roProfile->numRays(); ++rayIdx)
    {
      const Ray_& ray = roProfile->ray(rayIdx);

      //  Insert the lat-lon values, as well as the characteristic time of the profile
      //  for each sampled location associated with the profile.
      if (ray->assumeSphericalSymmetry()) {
        // All nodes share one location and time.
        lats.insert(lats.end(), ray->tpLat());
        lons.insert(lons.end(), ray->tpLon());
        times.insert(times.end(), extItr->time());
      } else {
        // Each node has a separate location, but the same time.
        lats.insert(lats.end(), ray->lats().cbegin(), ray->lats().cend());
        lons.insert(lons.end(), ray->lons().cbegin(), ray->lons().cend());
        times.insert(times.end(), ray->lats().size(), extItr->time());
      }
    }
    ++extItr;
  }

  if (lats.size() != times.size())
  {
    oops::Log::error() << funcName << "Number of locations and times are inconsistent"
            << ": lats (" << lats.size() << "), times (" << times.size()
            << "); they must be the same" << std::endl;
    assert(lats.size() == times.size());
  }

  if (lats.size() != numLocs)
  {
    oops::Log::error() << funcName << "Number of locations and computed numLocs are "
            << "inconsistent: lats (" << lats.size() << "), numLocs (" << numLocs
            << "); they must be the same" << std::endl;
    assert(lats.size() == times.size());
  }

  return numLocs;
}
// -----------------------------------------------------------------------------

void GnssroRayPathOrchestrator::fillPathsGroupedByLocation(
        std::size_t numObs,
        std::vector<util::Range<std::size_t>> & pathsGroupedByLocation) const
{
  const char funcName[]  = "GnssroRayPathOrchestrator::fillPathsGroupedByLocation";
  pathsGroupedByLocation.resize(numObs);

  std::size_t obIdx = 0;
  std::size_t locIdx = 0;
  for (Profiles_::const_iterator itr = profiles_.cbegin(); itr != profiles_.cend(); ++itr)
  {
    const Profile_& roProfile = *itr;
    for (std::size_t rayIdx = 0; rayIdx < roProfile->numRays(); ++rayIdx)
    {
      const Ray_& ray = roProfile->ray(rayIdx);

      // Do not go out of bounds of the array, but keep a count of the
      // number of rays in all the profiles.
      if (obIdx >= numObs)
      {
        ++obIdx;
        continue;
      }

      //  Set the range of sampled location indices used by this ray.
      //  It will have one sampled location for each lat-lon location in the ray.
      //  This is either 1 location per node, or a single node for the entire
      //  ray if the spherical symmetry assumption is in place.
      std::uint32_t rayNumLocs = ray->numLocations();
      pathsGroupedByLocation[obIdx].begin = locIdx;
      pathsGroupedByLocation[obIdx].end = locIdx + rayNumLocs;
      oops::Log::debug() << funcName << ": obIdx=" << obIdx << " of "
              << numObs << ", add " << rayNumLocs << " sampled locs "
              << "[" << pathsGroupedByLocation[obIdx].begin
              << "," << pathsGroupedByLocation[obIdx].end << ")" << std::endl;
      locIdx += rayNumLocs;
      ++obIdx;
    }
  }

  if (obIdx != numObs)
  {
    oops::Log::error() << funcName << ": found " << obIdx
            << " ray paths to specify in pathsGroupedByLocation, but there are "
            << numObs << " observations to specify; these must be the same."
            << std::endl;
  }
  return;
}
// -----------------------------------------------------------------------------

std::size_t GnssroRayPathOrchestrator::totalNumRays() const
{
  std::size_t numRays = 0;
  for (Profiles_::const_iterator p = profiles_.cbegin(); p != profiles_.cend(); ++p)
  {
    const Profile_& profile = *p;
    numRays += profile->numRays();
  }
  return numRays;
}
// -----------------------------------------------------------------------------

double GnssroRayPathOrchestrator::pressureInterpApply(
    int nlevs,
    const double * pres,
    const double * temp,
    int wi,
    double wf)
{
  // Interpolate pressure in a way that ensures both the bottom
  // and top of the layer are used. By contrast, the hydrostatic formulation
  // used in RefNCEP only uses pressure at the bottom level. Using pressure
  // at both bottom and top ensures better agreement between non-linear and
  // linear models if pressure is only perturbed at one of those levels.
  //
  // Ideally, we use a form of the hypsometric equation with a
  // non-zero lapse rate, but we fall back to isothermal interpolation
  // if temperature is effectively constant over the layer.
  // This formulation is similar to pressure interpolation used in ROPP.
  const double NEAR_ZERO = 1e-10;
  int cwi = wi - 1;  // convert 1-based index wi to 0-based cwi
  double interpPres = util::missingValue<double>();
  if (std::abs(temp[wi] - temp[cwi]) < NEAR_ZERO) {
    // layer is isothermal: use isothermal interpolation
    if (pres[cwi] <= 0.0 || pres[wi] <= 0.0) {
      // Min pressure <= 0.0: use linear interpolation
      interpPres = wf * pres[cwi] + (1.0 - wf) * pres[wi];
    } else {
      // Pressure is reasonable: use exponential interpolation
      interpPres = pres[cwi] * std::exp(wf * std::log(pres[wi]/pres[cwi]));
    }
  } else {
    // Use adapted hypsometric equation, scaled to temperature and pressure
    // values at the layer's bottom and top.
    double gamma = std::log(temp[wi]/temp[cwi])
                 / std::log(pres[wi]/pres[cwi]);
    double interpTemp = wf * temp[cwi] + (1.0 - wf) * temp[wi];
    interpPres = pres[cwi] * std::pow(interpTemp / temp[cwi], 1.0 / gamma);
  }
  return interpPres;
}

// -----------------------------------------------------------------------------
// Class method for computing refractivity from interpolated T, P, and Q.
//
double GnssroRayPathOrchestrator::modelRefractivity(
        float geomHgt, float lat, const GnssroGeoVaLs & ggv, std::size_t profileIdx,
        const std::unique_ptr<ufo::RefractivityCalculator> & refrCalc)
{
  // Get references to model vertical profiles (columns) for this profileIdx.
  const GnssroGeoVaLs::Profile_& gphProf = ggv.gphProfile(profileIdx);
  const GnssroGeoVaLs::Profile_& tempProf = ggv.tempProfile(profileIdx);
  const GnssroGeoVaLs::Profile_& sphumProf = ggv.sphumProfile(profileIdx);
  const GnssroGeoVaLs::Profile_& presProf = ggv.presProfile(profileIdx);

  // Convert geometric height to geopotential height
  float gph = formulas::Geometric_to_Geopotential_Height(lat, geomHgt);

  // Compute vertical interpolation weights within a model profile for the
  // target geopotential height
  int wi;     // Base vertical index for vertical interpolation (1-based).
  double wf;  // Weight for base index (rest of weight goes to idx wi + 1).
  vert_interp_weights_f90(ggv.nlevs(), gph, gphProf.data(), wi, wf);

  // Linearly interpolate temp and sphum to target gph
  double temp;    // Interpolated air temperature in K
  double sphum;   // Interpolated specific humidity in kg/kg
  vert_interp_apply_f90(ggv.nlevs(), tempProf.data(), temp, wi, wf);
  vert_interp_apply_f90(ggv.nlevs(), sphumProf.data(), sphum, wi, wf);
  double pres = pressureInterpApply(
          ggv.nlevs(), presProf.data(), tempProf.data(), wi, wf);

  // Compute model refractivity from the interpolated temp,sphum,pres.
  double refr = refrCalc->N(temp, sphum, pres);
  return refr;
}

// -----------------------------------------------------------------------------
// Class method for computing derivatives of refractivity from interpolated T, P, and Q.
GnssroRayPathOrchestrator::TrajTuple
GnssroRayPathOrchestrator::modelTrajectory(
        float geomHgt, float lat, const GnssroGeoVaLs & ggv, std::size_t profileIdx,
        const std::unique_ptr<ufo::RefractivityCalculator> & refrCalc)
{
  // Get references to model vertical profiles (columns) for this profileIdx.
  const GnssroGeoVaLs::Profile_& gphProf = ggv.gphProfile(profileIdx);
  const GnssroGeoVaLs::Profile_& tempProf = ggv.tempProfile(profileIdx);
  const GnssroGeoVaLs::Profile_& sphumProf = ggv.sphumProfile(profileIdx);
  const GnssroGeoVaLs::Profile_& presProf = ggv.presProfile(profileIdx);

  // Convert geometric height to geopotential height
  float gph = formulas::Geometric_to_Geopotential_Height(lat, geomHgt);

  // Compute vertical interpolation weights within a model profile for the
  // target geopotential height
  int wi;     // Base vertical index for vertical interpolation (1-based).
  double wf;  // Weight for base index (rest of weight goes to idx wi + 1).
  vert_interp_weights_f90(ggv.nlevs(), gph, gphProf.data(), wi, wf);

  // Linearly interpolate temp and sphum to target gph
  double temp;    // Interpolated air temperature in K
  double sphum;   // Interpolated specific humidity in kg/kg
  vert_interp_apply_f90(ggv.nlevs(), tempProf.data(), temp, wi, wf);
  vert_interp_apply_f90(ggv.nlevs(), sphumProf.data(), sphum, wi, wf);
  double pres = pressureInterpApply(
          ggv.nlevs(), presProf.data(), tempProf.data(), wi, wf);

  // Compute derivatives of model refractivity from the interpolated temp,sphum,pres.
  double dNdT = refrCalc->dNdT(temp, sphum, pres);
  double dNdP = refrCalc->dNdP(temp, sphum, pres);
  double dNdQ = refrCalc->dNdQ(temp, sphum, pres);

  // Return the trajectory elements as a tuple.
  TrajTuple traj(dNdT, dNdQ, dNdP, wf, wi);
  return traj;
}

// -----------------------------------------------------------------------------
// Associated ostream operator
std::ostream& operator<<(std::ostream& os, const GnssroRayPathOrchestrator& orc)
{
  //  Accumulate stats over all the profiles managed by the orchestrator
  const GnssroRayPathOrchestrator::Profiles_& profiles = orc.profiles();
  std::size_t numProfiles = profiles.size();  // Total RO profiles/seqNums
  std::size_t numRays = 0;         // Total observation samples turned into rays.
  std::size_t numNodes = 0;
  std::size_t numSampledLocations = 0;
  std::size_t maxNodesForRay = 0;
  int maxNodeSeqNum = 0;
  std::size_t maxNodeRayIdx = 0;

  for (GnssroRayPathOrchestrator::Profiles_::const_iterator itr = profiles.cbegin();
       itr != profiles.cend(); ++itr)
  {
    const GnssroRayPathOrchestrator::Profile_& roProfile = *itr;
    numRays += roProfile->numRays();
    for (std::size_t rayIdx = 0; rayIdx < roProfile->numRays(); ++rayIdx)
    {
      const GnssroRayPathOrchestrator::Ray_& ray = roProfile->ray(rayIdx);

      numSampledLocations += ray->numLocations();
      numNodes += ray->numNodes();
      if (ray->numNodes() > maxNodesForRay) {
        maxNodesForRay = ray->numNodes();
        maxNodeSeqNum = roProfile->seqNum();
        maxNodeRayIdx = rayIdx;
      }
    }
  }

  // Ratio of sampledLocations (columns from the model data) to the number of
  // tangent point observations. A 1D operator would return 1.0, so this is
  // a metric of how efficiently we are sampling the model data.
  float sampledLocationsPerRay = static_cast<double>(numSampledLocations)
                               / static_cast<double>(numRays);

  os << "ROprofiles=" << numProfiles
     << ", rays=" << numRays
     << ", numSampLocs=" << numSampledLocations
     << ", sampLocs/ray=" << sampledLocationsPerRay
     << ", numRayNodes=" << numNodes
     << ", maxNodesForRay=" << maxNodesForRay
     << " (seqNum=" << maxNodeSeqNum
     << ", rayIdx=" << maxNodeRayIdx << ")";
  return os;
}
// -----------------------------------------------------------------------------

}  // namespace ufo

