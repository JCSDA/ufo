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
#include "ufo/filters/VariableTransformParametersBase.h"
#include "ufo/operators/gnssro/utils/RefractivityCalculator.h"
#include "ufo/utils/Constants.h"
#include "ufo/utils/VertInterp.interface.h"
#include "ufo/variabletransforms/Cal_NonLocalPseudoExcessPhase.h"

namespace ufo {

/************************************************************************************/
//  Cal_NonLocalPseudoExcessPhase
/************************************************************************************/

// Register the variable transform with ufo::TransformFactory
static TransformMaker<Cal_NonLocalPseudoExcessPhase>
    makerCal_NonLocalPseudoExcessPhase("NonLocalPseudoExcessPhase");

Cal_NonLocalPseudoExcessPhase::Cal_NonLocalPseudoExcessPhase(
    const Parameters_ &options,
    const ObsFilterData &data,
    const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
    const std::shared_ptr<ioda::ObsDataVector<float>> &obserr)
  : TransformBase(options, data, flags, obserr)
  , group_(options.group)
  , refrvariable_(options.refrVariable)
  , nlpepvariable_(options.nlpepVariable)
  , rpParams_(options.rayPathGenType, options.approxRayLength,
              options.top2D, options.res, options.nHoriz)
  , orchestrator_(obsdb_, rpParams_)
{}

/************************************************************************************/

void Cal_NonLocalPseudoExcessPhase::runTransform(const std::vector<bool> &apply) {
  const char funcName[] = "Cal_NonLocalPseudoExcessPhase::runTransform";
  oops::Log::trace() << funcName << ": START: obsName: " << obsName()
                     << ", numProfiles=" << orchestrator_.profiles().size()
                     << ", nlocs=" << obsdb_.nlocs() << std::endl;

  // Confirm our sizing for number of obs is what we expect.
  const size_t nlocs = obsdb_.nlocs();
  std::size_t totalNumRays = orchestrator_.totalNumRays();
  if (totalNumRays != nlocs)
  {
    oops::Log::error() << funcName << ": size mismatch between obsdb nlocs ("
            << nlocs << ") and Orchestrator's total number of rays ("
            << totalNumRays << ")" << std::endl;
    throw eckit::BadValue("size mismatch between obsdb nlocs and Orchestrator in "
                          "Cal_NonLocalPseudoExcessPhase", Here());
  }

  // Allocate space for nonLocalPseudoExcessPhase
  std::vector<float> nlpep(nlocs, missingValueFloat);

  // Get observed refractivity
  std::vector<double> refractivities;
  std::vector<double> geomHeights;
  getObservation(group_, refrvariable_, refractivities, true);
  getObservation("MetaData", "height", geomHeights, true);

  double dmiss = util::missingValue<double>();

  std::size_t ovecIdx = 0;  // Index into locations in ObsVector (nlocs)
  //  Iterate over RO profiles.
  const GnssroRayPathOrchestrator::Profiles_& profiles = orchestrator_.profiles();
  for (const auto& roProfile : profiles)
  {
    int seqNum = roProfile->seqNum();
    oops::Log::debug() << funcName << ": processing profile seqNum=" << seqNum
            << " with " << roProfile->numRays() << " rays" << std::endl;

    // Get subset of the observation arrays associated with this RO profile.
    const double * profileGeomHeights = geomHeights.data() + ovecIdx;
    const double * profileRefractivities = refractivities.data() + ovecIdx;

    // Iterate over rays within this profile.
    for (std::size_t rayIdx = 0; rayIdx < roProfile->numRays(); ++rayIdx, ++ovecIdx)
    {
      const GnssroRayPathOrchestrator::Ray_& ray = roProfile->ray(rayIdx);

      // Do not compute NLPEP if this obs has been excluded by the where statement
      if (!apply[ovecIdx]) continue;

      // If the tangent point refractivity is missing, skip this ob.
      if (refractivities[ovecIdx] == dmiss) continue;

      std::size_t tpNodeIdx = ray->tpNodeIdx();
      double accumNLPEP = 0.0;
      for (std::size_t nodeIdx = 0; nodeIdx < ray->numNodes(); ++nodeIdx)
      {
        // Get refractivity at the average node height, assuming a
        // spherically symmetric atmosphere. This allows us to treat
        // the observations in the profile as though they formed a
        // vertical column.
        double refr = dmiss;
        double geomHgt = dmiss;
        if (nodeIdx == tpNodeIdx) {
          // Tangent point of ray: use the observed values.
          geomHgt = geomHeights[ovecIdx];
          refr = refractivities[ovecIdx];
        } else {   // Average of heights of edges of this ray segment.
          float thisGeomHgt = ray->geomHgt(nodeIdx);
          float nextGeomHgt = ray->geomHgt(nodeIdx + 1);
          geomHgt = static_cast<double>(0.5 * (thisGeomHgt + nextGeomHgt));

          // Get vertical interpolation weights using geometric height.
          // Implementation note: If the height is out-of-bounds, this will give us
          // weights for 100% of the topmost or bottommost value.
          int wi;     // Base vertical index for vertical interpolation (1-based).
          double wf;  // Weight for base index (rest of weight goes to idx wi + 1).
          vert_interp_weights_f90(roProfile->numRays(), geomHgt, profileGeomHeights, wi, wf);

          // Vertically interpolate refractivity to target geometric height
          vert_interp_apply_f90(roProfile->numRays(), profileRefractivities, refr, wi, wf);
        }

        double segmentLength = static_cast<double>(ray->segLen(nodeIdx));
        if (refr != dmiss) {
          accumNLPEP += (refr * segmentLength);
        } else {
          oops::Log::warning() << funcName << ": seqNum=" << seqNum << ", rayIdx=" << rayIdx
                  << ", ovecIdx=" << ovecIdx << ", nodeIdx=" << nodeIdx << ", tpNodeIdx="
                  << tpNodeIdx << ", numNodes=" << ray->numNodes() << ", geomHgt=" << geomHgt
                  << ", tpRefr=" << refractivities[ovecIdx]
                  << ", refr is missing for this node, so NLPEP set to missing"
                  << std::endl;
          accumNLPEP = dmiss;
          break;
        }
      }  // end of iteration over nodes in a single ray

      if (accumNLPEP != dmiss) {
        // All the nodes had non-missing data, so we could get a complete accumulation.
        // The factor converts units from Refractivity (N) to excess index of refraction.
        // The resulting units are in meters.
        nlpep[ovecIdx] = accumNLPEP * RefractivityCalculator::N_TO_EXCESS_IOR;
      }
    }  // end iteration over rays in a single RO profile
  }  // end iteration over RO profiles

  // put new variable at existing locations
  putObservation(nlpepvariable_, nlpep, getDerivedGroup(group_));
  oops::Log::trace() << funcName << ": DONE: obsName: " << obsName() << ", added "
                     << nlpep.size() << " obs to " << nlpepvariable_ << " in group "
                     << getDerivedGroup(group_) << std::endl;
  return;
}

}  // namespace ufo
