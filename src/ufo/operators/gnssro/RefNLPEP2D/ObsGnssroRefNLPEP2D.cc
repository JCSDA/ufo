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
#include <ostream>
#include <utility>
#include <vector>

#include "eckit/exception/Exceptions.h"
#include "ioda/ObsVector.h"

#include "oops/base/Locations.h"
#include "oops/base/Variables.h"
#include "oops/interface/SampledLocations.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/ObsTraits.h"
#include "ufo/operators/gnssro/RefNLPEP2D/ObsGnssroRefNLPEP2D.h"
#include "ufo/operators/gnssro/utils/GnssroGeoVaLs.h"
#include "ufo/operators/gnssro/utils/GnssroProfileExtractor.h"
#include "ufo/operators/gnssro/utils/GnssroRayPathGenerator.h"
#include "ufo/operators/gnssro/utils/RefractivityCalculator.h"
#include "ufo/SampledLocations.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsGnssroRefNLPEP2D> makerGnssroRefNLPEP2D_(
    "GnssroRefNLPEP2D");
// -----------------------------------------------------------------------------

ObsGnssroRefNLPEP2D::ObsGnssroRefNLPEP2D(const ioda::ObsSpace & odb,
                       const Parameters_ & parameters)
  : ObsOperatorBase(odb)
  , odb_(odb)
  , varin_()
  , opts_(parameters.options.value())
  , rpParams_(opts_.rayPathGenType, opts_.approxRayLength,
              opts_.top2D, opts_.res, opts_.nHoriz)
  , orchestrator_(odb_, rpParams_)
  , refrCalc_(RefractivityCalculator::create(opts_.refrAlgorithm, opts_.useCompress))
{
  oops::Log::trace() << "ObsGnssroRefNLPEP2D - c'tor body starting" << std::endl;
  const std::vector<std::string> vv{GnssroGeoVaLs::VAR_TEMP, GnssroGeoVaLs::VAR_SPHUM,
          GnssroGeoVaLs::VAR_PRES, GnssroGeoVaLs::VAR_GPH};
  varin_.reset(new oops::Variables(vv));

  // Provide aggregated information about the ray creation process
  oops::Log::info() << "ObsGnssroRefNLPEP2D ctor summary: " << orchestrator_ << std::endl;

  oops::Log::trace() << "ObsGnssroRefNLPEP2DParameters: " << parameters.toConfiguration()
                       << std::endl;
  oops::Log::trace() << "ObsGnssroRefNLPEP2D created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsGnssroRefNLPEP2D::~ObsGnssroRefNLPEP2D() {
  oops::Log::trace() << "ObsGnssroRefNLPEP2D destroyed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGnssroRefNLPEP2D::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                             ObsDiagnostics & d, const QCFlags_t & qc_flags) const {
  const char funcName[] = "ObsGnssroRefNLPEP2D::simulateObs";
  oops::Log::trace() << funcName << " started" << std::endl;

  oops::Log::debug() << funcName << ": GeoVaLs has " << gv.nlocs() << " locations and "
          << gv.getVars().size() << " variables:" << std::endl;
  for (std::size_t varIdx = 0; varIdx < gv.getVars().size(); ++varIdx)
  {
    const oops::Variable& var = (gv.getVars())[varIdx];
    oops::Log::debug() << var.name() << ": " << gv.nlevs(var) << " levels, "
          << gv.nprofiles(var) << " profiles/paths" << std::endl;
  }
  oops::Log::debug() << funcName << ": ovec has " << ovec.nlocs() << " locations, "
      << ovec.size() << " array elements, " << ovec.nvars() << " variables" << std::endl;

  // Confirm our sizing for number of obs is what we expect.
  std::size_t totalNumRays = orchestrator_.totalNumRays();
  if (totalNumRays != ovec.nlocs())
  {
    oops::Log::error() << funcName << ": size mismatch between ObsVector ("
            << ovec.nlocs() << ") and Orchestrator's total number of rays ("
            << totalNumRays << ")" << std::endl;
    throw eckit::BadValue("size mismatch between ObsVector and Orchestrator in "
                          "ObsGnssroRefNLPEP2D", Here());
  }

  // Return early if there are no geoVals to operate on.
  if (gv.nlocs() == 0)
  {
    return;
  }

  // Extract the model data from the GeoVaLs.
  GnssroGeoVaLs ggv(gv);

  // Get geometric height and latitude from the observations.
  std::vector<float> latitudes(odb_.nlocs());
  std::vector<float> longitudes(odb_.nlocs());
  std::vector<double> geomHeights(odb_.nlocs());
  odb_.get_db("MetaData", "latitude", latitudes);
  odb_.get_db("MetaData", "longitude", longitudes);
  odb_.get_db("MetaData", "height", geomHeights);

  //  Iterate over RO profiles.
  double dmiss = util::missingValue<double>();
  std::size_t ovecIdx = 0;
  std::size_t profileIdx = 0;
  const GnssroRayPathOrchestrator::Profiles_& profiles = orchestrator_.profiles();
  for (const auto& roProfile : profiles)
  {
    int seqNum = roProfile->seqNum();

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
      auto fwdItr = rayItr.base();
      std::size_t rayIdx = std::distance(roProfile->cbegin(), fwdItr) - 1;
      const GnssroProfileRayPaths::Ray_& ray = *rayItr;

      // Get the model profile range for this ray
      std::size_t rayEndProfileIdx = profileIdx;
      std::size_t rayStartProfileIdx = rayEndProfileIdx - ray->numLocations();
      profileIdx = rayStartProfileIdx;
      oops::Log::debug() << funcName << ": seqNum=" << seqNum << ", rayIdx= "
          << rayIdx << ", numSampLocs=" << ray->numLocations()
          << ", ray profileIdx range=[" << rayStartProfileIdx
          << "," << rayEndProfileIdx << "), starting at profileIdx="
          << profileIdx << std::endl;

      std::size_t tpNodeIdx = ray->tpNodeIdx();

      oops::Log::debug() << funcName << ": seqNum=" << seqNum << ", rayIdx=" << rayIdx
             << ", ovecIdx=" << ovecIdx << ", tpNodeIdx=" << tpNodeIdx
             << ", computing NLPEP over numNodes=" << ray->numNodes()
             << ", totalLength=" << ray->totalLength()
             << ", tpHgt=" << geomHeights[ovecIdx] << std::endl;

      //  Iterate over nodes in the ray. We accumulate:
      //    fwdNLPEP = Model 2D Non-local pseudo excess phase. This is the forward part
      //               from integrating N(x,y) over the ray path.
      double fwdNLPEP = 0.0;
      for (std::size_t nodeIdx = 0; nodeIdx < ray->numNodes(); ++nodeIdx)
      {
        // Unnormalized segment length for this node of the ray in meters
        double segmentLength = ray->segLen(nodeIdx);

        // Get the geometric height of the node and observed latitude of tangent point.
        // Compute corresponding geopotential height of this observation sample.
        //
        // Use the average of the edge heights at all locations except the tangent point.
        // The edge average will not find the lowest height at the tangent point.
        // Instead, it will use two edges that are both higher than the tangent point.
        // So we use the observed geometric height for the tangent point instead.
        //
        double avgGeomHgt;
        if (nodeIdx == tpNodeIdx) {
          avgGeomHgt = geomHeights[ovecIdx];
        } else {
          float thisGeomHgt = ray->geomHgt(nodeIdx);
          float nextGeomHgt = ray->geomHgt(nodeIdx + 1);
          avgGeomHgt = static_cast<double>(0.5 * (thisGeomHgt + nextGeomHgt));
        }
        double refr = GnssroRayPathOrchestrator::modelRefractivity(
                avgGeomHgt, latitudes[ovecIdx], ggv, profileIdx, refrCalc_);

        // Integrate refractivity by the segment length of this ray node
        // to get nonLocalPseudoExcessPhase contribution from the horizontally
        // varying model data.
        double fwdNLPEPincr = (refr * segmentLength);
        fwdNLPEP += fwdNLPEPincr;
        oops::Log::debug() << funcName << ": seqNum=" << seqNum << ", rayIdx=" << rayIdx
              << ", nodeIdx=" << nodeIdx << ", profileIdx=" << profileIdx
              << ", geomHgt=" << avgGeomHgt
              << ", increment fwdNLPEP=" << fwdNLPEP << ", refr=" << refr
              << ", segmentLength=" << segmentLength
              << ", NLPEPincr=" << fwdNLPEPincr << std::endl;

        oops::Log::debug() << funcName << ": seqNum=" << seqNum << ", rayIdx=" << rayIdx
             << ", nodeIdx=" << nodeIdx << ", ovecIdx=" << ovecIdx
             << ": Node complete: refr=" << refr << ", segLen=" << segmentLength
             << ", fwdNLPEP=" << fwdNLPEP << std::endl;

        // Only change location for next node if we are not assuming
        // spherical symmetry, which means all nodes have the same location.
        if (!ray->assumeSphericalSymmetry()) {
          ++profileIdx;
        }
      }  // end iteration over nodes in a ray

      // Scale from refractivity N to excess index of refraction (n-1).
      // This is a unit conversion from N * m to m.
      ovec[ovecIdx] = fwdNLPEP * RefractivityCalculator::N_TO_EXCESS_IOR;

      // Prepare for the next lower ray by setting the profileIdx to point to
      // the start of the model profiles for the current ray, which is also the
      // end of the range for the next lower ray within the profile.
      profileIdx = rayStartProfileIdx;
    }   // end iteration over rays in a single RO profile

    // Update ob vector index to skip to the next profile.
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
  }  // end iteration over RO profiles.

  oops::Log::trace() << funcName << " complete" << std::endl;
}
// -----------------------------------------------------------------------------
ObsGnssroRefNLPEP2D::Locations_
ObsGnssroRefNLPEP2D::locations() const {
    typedef oops::SampledLocations<ObsTraits> SampledLocations_;
    oops::Log::trace() << "ObsGnssroRefNLPEP2D::locations() enter" << std::endl;

  // Get the arrays of latitudes, longitudes, and datetimes used to create the
  // sampled locations. This method will allocate the vectors properly.
  std::vector<float> lons;
  std::vector<float> lats;
  std::vector<util::DateTime> times;
  std::size_t numSampledLocs = orchestrator_.getSampledLatLonTimes(lats, lons, times);

  // Set range of sampled locations used by each observed ray in each profile.
  std::vector<util::Range<std::size_t>> pathsGroupedByLocation;
  orchestrator_.fillPathsGroupedByLocation(odb_.nlocs(), pathsGroupedByLocation);

  SampledLocations_ sampledLocations(
          std::make_unique<SampledLocations>(lons, lats, times, odb_.distribution(),
                                             std::move(pathsGroupedByLocation)));
  oops::Log::debug() << "ObsGnssroRefNLPEP2D::locations: creating SampledLocations for "
                     << odb_.nlocs() << " ob samples to make " << numSampledLocs
                     << " model paths; SampledLocs/ob="
                     << static_cast<float>(numSampledLocs)/static_cast<float>(odb_.nlocs())
                     << ", locationsSampledOnceAndInOrder="
                     << sampledLocations.sampledLocations().areLocationsSampledOnceAndInOrder()
                     << std::endl;

  const std::vector<util::Range<std::size_t>>& pathRanges =
          sampledLocations.sampledLocations().pathsGroupedByLocation();
  if (pathRanges.size() != odb_.nlocs())
  {
    oops::Log::error() << "ObsGnssroRefNLPEP2D::locations: SampledLocations has "
        << pathRanges.size() << " entries, but must have one range for each of "
        << odb_.nlocs() << " observations" << std::endl;
  }

  std::vector<int> seqNums(odb_.nlocs());
  odb_.get_db("MetaData", "sequenceNumber", seqNums);

  for (std::size_t obIdx = 0; obIdx < odb_.nlocs(); ++obIdx)
  {
    const util::Range<size_t> & pathRange = pathRanges[obIdx];
    oops::Log::debug() << "SampledLoc for obIdx " << obIdx << " in profile seqNum="
        << seqNums[obIdx] << " has " << pathRange.end - pathRange.begin
        << " locations from " << pathRange.begin << " to " << pathRange.end << "; "
        << "lat/lon[" << pathRange.begin << "]=("
        << (sampledLocations.latitudes())[pathRange.begin]
        << "," << (sampledLocations.longitudes())[pathRange.begin] << "), "
        << "lat/lon[" << pathRange.end - 1 << "]=("
        << (sampledLocations.latitudes())[pathRange.end-1]
        << "," << (sampledLocations.longitudes())[pathRange.end-1] << ")" << std::endl;
  }
  return sampledLocations;
}

// -----------------------------------------------------------------------------

void ObsGnssroRefNLPEP2D::computeReducedVars(const oops::Variables & reducedVars,
                                            GeoVaLs & /*geovals*/) const {
  // No method for reducing the set of nhoriz_ profiles sampling each location into a single
  // profile has been implemented so far, so when this obs operator is in use, neither it nor
  // any obs filters or bias predictors can request variables in the reduced format.
  if (reducedVars.size() != 0)
    throw eckit::NotImplemented("ObsGnssroRefNLPEP2D is unable to compute the reduced "
                                "representation of GeoVaLs", Here());
}

// -----------------------------------------------------------------------------

void ObsGnssroRefNLPEP2D::print(std::ostream & os) const {
  os << "ObsGnssroRefNLPEP2D::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
