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
#include <algorithm>
#include "ufo/filters/VariableTransformParametersBase.h"
#include "ufo/operators/gnssro/utils/RefractivityCalculator.h"
#include "ufo/utils/Constants.h"
#include "ufo/utils/VertInterp.interface.h"
#include "ufo/variabletransforms/Cal_GnssroRefractivityGradient.h"

namespace ufo {

/************************************************************************************/
//  Cal_GnssroRefractivityGradient
/************************************************************************************/

// Register the variable transform with ufo::TransformFactory
static TransformMaker<Cal_GnssroRefractivityGradient>
    makerCal_GnssroRefractivityGradient("GnssroRefractivityGradient");

Cal_GnssroRefractivityGradient::Cal_GnssroRefractivityGradient(
    const Parameters_ &options,
    const ObsFilterData &data,
    const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
    const std::shared_ptr<ioda::ObsDataVector<float>> &obserr)
  : TransformBase(options, data, flags, obserr)
  , group_(options.group)
  , refrvariable_(options.refrVariable)
  , heightvariable_(options.heightVariable)
  , gradvariable_(options.gradVariable)
  , minSuperRefraction_(options.minSuperRefraction)
  , minSRThreshold_(minSuperRefraction_ * Constants::superRefractionCritVal)  // N/m
  , nearThreshold_(0.6 * Constants::superRefractionCritVal)  // N/m
  , quiteNearThreshold_(0.7 * Constants::superRefractionCritVal)  // N/m
  , veryNearThreshold_(0.8 * Constants::superRefractionCritVal)  // N/m
  , extremelyNearThreshold_(0.9 * Constants::superRefractionCritVal)  // N/m
  , calcDuctingFlag_(options.calculateDuctingFlag)
  , calcProfileDuctingFlag_(options.calculateProfileDuctingFlag)
  , extractor_(obsdb_)
{
  // Validate configured setting of min super refraction
  if (minSuperRefraction_ < 0.0 || minSuperRefraction_ >= 0.6) {
    throw eckit::BadValue("GnssroRefractivityGradient: bad min super refraction value "
      + std::to_string(minSuperRefraction_)
      + ": must be >= 0 and < near ducting threshold of 0.6");
  }
}

/************************************************************************************/

void Cal_GnssroRefractivityGradient::runTransform(const std::vector<bool> &apply) {
  const char funcName[] = "Cal_GnssroRefractivityGradient::runTransform";
  oops::Log::trace() << funcName << ": START: obsName: " << obsName() << std::endl;

  // Allocate space for refractivity gradient and optional flag variables.
  const size_t nlocs = obsdb_.nlocs();
  std::vector<float> refrGradients(nlocs, missingValueFloat);

  std::vector<int> ductingFlags;
  if (calcDuctingFlag_) {
    ductingFlags.assign(nlocs, NONDUCTING);
  }
  std::vector<int> profileDuctingFlags;
  if (calcProfileDuctingFlag_) {
    profileDuctingFlags.assign(nlocs, NONDUCTING);
  }

  // Get observed refractivity and heights.
  std::vector<double> refractivities;
  std::vector<double> geomHeights;
  getObservation(group_, refrvariable_, refractivities, true);
  getObservation("MetaData", heightvariable_, geomHeights, true);

  double dmiss = util::missingValue<double>();

  // Iterate over GNSSRO profiles found by the GnssroProfileExtractor.
  for (GnssroProfileExtractor::const_iterator itr = extractor_.cbegin();
       itr != extractor_.cend(); ++itr)
  {
    oops::Log::debug() << funcName << ": Profile slice "
                       << itr->seqNum() << " from " << itr->start()
                       << " to " << itr->end() << std::endl;

    // Iterate over all the levels in the profile, excluding the topmost.
    // dN/dz for a given level will be determined by the differences
    // between that level and the level above it.
    double maxDuctingFlag = NONDUCTING;
    std::size_t topObIdx = itr->end() - 1;
    for (std::size_t obIdx = itr->start(); obIdx < topObIdx; ++obIdx)
    {
      // Do not compute derived values if ob has been excluded by the where statement
      if (!apply[obIdx]) continue;

      // If the refractivity or height is missing, skip this ob.
      if (refractivities[obIdx] == dmiss) continue;
      if (geomHeights[obIdx] == dmiss) continue;

      // Compute dN/dZ in units N-units / meter
      refrGradients[obIdx] = (refractivities[obIdx+1] - refractivities[obIdx])
                           / (geomHeights[obIdx+1] - geomHeights[obIdx]);

      oops::Log::debug() << funcName << ", seqNum=" << itr->seqNum()
         << ", obIdx=" << obIdx << ", Ntop=" << refractivities[obIdx+1]
         << ", Nbot=" << refractivities[obIdx]
         << ", Htop=" << geomHeights[obIdx+1]
         << ", Hbot=" << geomHeights[obIdx]
         << ", dN/dz=" << refrGradients[obIdx]
         << ", minSRThreshold=" << minSRThreshold_ << std::endl;

      // Set optional flags indicating non-ducting, near-ducting, or ducting
      if (calcDuctingFlag_ || calcProfileDuctingFlag_) {
        double absRefGradient = std::abs(refrGradients[obIdx]);
        int ductingFlag = NONDUCTING;
        if (absRefGradient >= Constants::superRefractionCritVal) {
          ductingFlag = DUCTING;
        } else if (absRefGradient >= extremelyNearThreshold_) {
          ductingFlag = EXTREMELY_NEAR_DUCTING;
        } else if (absRefGradient >= veryNearThreshold_) {
          ductingFlag = VERY_NEAR_DUCTING;
        } else if (absRefGradient >= quiteNearThreshold_) {
          ductingFlag = QUITE_NEAR_DUCTING;
        } else if (absRefGradient >= nearThreshold_) {
          ductingFlag = NEAR_DUCTING;
        } else if (absRefGradient >= minSRThreshold_) {
          ductingFlag = MIN_SUPER_REFRACTION;
        }

        if (ductingFlag != 0) {
          oops::Log::debug() << funcName << ", seqNum=" << itr->seqNum()
              << ", obIdx=" << obIdx << ", detected ductingFlag="
              << ductingFlag << std::endl;
        }

        if (calcDuctingFlag_) {
          ductingFlags[obIdx] = ductingFlag;
        }
        if (calcProfileDuctingFlag_) {
          if (ductingFlag > maxDuctingFlag) {
            maxDuctingFlag = ductingFlag;
          }
        }
      }  // end logic for computing ducting flags.
    }  // end iteration over observations in a single RO profile

    // Set ducting flag for the entire profile based on the level with the
    // more extreme ducting properties.
    if (calcProfileDuctingFlag_ && maxDuctingFlag != NONDUCTING) {
        std::fill(profileDuctingFlags.begin() + itr->start(),
                  profileDuctingFlags.begin() + itr->end(),
                  maxDuctingFlag);
    }
    oops::Log::debug() << funcName << ", seqNum=" << itr->seqNum()
      << ", profileDuctingFlag=" << maxDuctingFlag << std::endl;
  }  // end iteration over RO profiles

  // put new variable at existing locations
  putObservation(gradvariable_, refrGradients, getDerivedGroup(group_));
  if (calcDuctingFlag_) {
    putObservation("sampleDuctingFlag", ductingFlags, getDerivedGroup(group_));
  }
  if (calcProfileDuctingFlag_) {
    putObservation("profileDuctingFlag", profileDuctingFlags, getDerivedGroup(group_));
  }
  oops::Log::trace() << funcName << ": DONE: obsName: " << obsName() << ", added "
                     << refrGradients.size() << " obs to " << gradvariable_ << " in group "
                     << getDerivedGroup(group_) << std::endl;
  return;
}

}  // namespace ufo
