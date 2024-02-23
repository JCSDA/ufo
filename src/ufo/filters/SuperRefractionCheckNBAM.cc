/*
 * (C) Copyright 2017-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/SuperRefractionCheckNBAM.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/QCflags.h"
#include "ufo/utils/Constants.h"
#include "ufo/variabletransforms/Formulas.h"

namespace ufo {

// -----------------------------------------------------------------------------

SuperRefractionCheckNBAM::SuperRefractionCheckNBAM(
                                 ioda::ObsSpace & obsdb,
                                 const Parameters_ & parameters,
                                 std::shared_ptr<ioda::ObsDataVector<int> > flags,
                                 std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, parameters, flags, obserr), parameters_(parameters)
{
  oops::Log::trace() << "SuperRefractionCheckNBAM"<< std::endl;
  allvars_ += Variable("ObsDiag/atmosphericRefractivity_model");
  allvars_ += Variable("ObsDiag/geopotentialHeight_model");
  allvars_ += Variable("MetaData/impactParameterRO");
  allvars_ += Variable("MetaData/earthRadiusCurvature");
  allvars_ += Variable("MetaData/latitude");
  allvars_ += Variable("MetaData/geoidUndulation");
  allvars_ += Variable("ObsValue/bendingAngle");
  // get the list of obs grouping variables, and check that sequenceNumber,
  // and only sequenceNumber is being used
  auto obsGroupVars = obsdb.obs_group_vars();
  if ((obsGroupVars.size() != 1) || (obsGroupVars[0] != "sequenceNumber")) {
    throw eckit::BadParameter("sequenceNumber and only sequenceNumber should be used"
    " in obs grouping", Here());
  }
  if (parameters_.step.value() < 1)
    throw eckit::BadParameter("step has to be equal to or larger than 1", Here());
  if (parameters_.checkHeight.value() < 2000.0)
    throw eckit::BadParameter("the check impact height has to be greater than 2km", Here());
  if (parameters_.closeLayers.value() < 1)
    throw eckit::BadParameter("the number of model layers close-to-SR "
                              "has to be equal to or larger than 1",  Here());
}
// -----------------------------------------------------------------------------

SuperRefractionCheckNBAM::~SuperRefractionCheckNBAM() {
  oops::Log::trace() << "SuperRefractionCheckNBAM: destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void SuperRefractionCheckNBAM::applyFilter(
                                      const std::vector<bool> & apply,
                                      const Variables & filtervars,
                                      std::vector<std::vector<bool>> & flagged) const {
  oops::Log::trace() << "SuperRefractionCheckNBAM postFilter " << std::endl;
  const float missingFloat = util::missingValue<float>();
  // Get the refractivity from the obs diagnostics, including the number of
  // vertical levels on which the refractivity has been calculated (nRefLevels)
  Variable refractivityVariable = Variable("ObsDiag/atmosphericRefractivity_model");
  oops::Log::debug() << data_.nlevs(refractivityVariable) << std::endl;
  const size_t nRefLevels = data_.nlevs(refractivityVariable);
  std::vector<std::vector<float>> refractivity;

  for (size_t iLevel = 0; iLevel < nRefLevels; ++iLevel) {
    std::vector<float> inputData;
    data_.get(refractivityVariable, static_cast<int>(iLevel), inputData);
    refractivity.push_back(inputData);
  }

  // Get the height of the levels on which the refractivity has been calculated.
  // Must be the same length as the array defining the refractivity.
  Variable modelHeightsVariable = Variable("ObsDiag/geopotentialHeight_model");
  oops::Log::debug() << data_.nlevs(modelHeightsVariable) << std::endl;
  if (data_.nlevs(modelHeightsVariable) != nRefLevels) {
    throw eckit::BadValue("Model heights and refractivity must have the same number of levels",
                          Here());
  }
  std::vector<std::vector<float>> modelHeights;

  // Read the heights of the refractivity levels
  for (size_t iLevel = 0; iLevel < nRefLevels; ++iLevel) {
    std::vector<float> inputData;
    data_.get(modelHeightsVariable, static_cast<int>(iLevel), inputData);
    modelHeights.push_back(inputData);
  }

  // For the benefits of debugging, output model height  for the first observation
  oops::Log::debug() << "model height of first obs";
  for (size_t iLevel = 0; iLevel < nRefLevels; ++iLevel) {
    oops::Log::debug() << modelHeights[iLevel][0]<< " ";
  }
  oops::Log::debug() << std::endl;

  ioda::ObsDataVector<float> impactParameterObs(obsdb_, "impactParameterRO", "MetaData");
  ioda::ObsDataVector<float> radiusCurvature(obsdb_, "earthRadiusCurvature", "MetaData");
  ioda::ObsDataVector<float> latitude(obsdb_, "latitude", "MetaData");
  ioda::ObsDataVector<float> geoid(obsdb_, "geoidUndulation", "MetaData");
  ioda::ObsDataVector<float> obsValue(obsdb_, "bendingAngle", "ObsValue");

  // Get the record numbers from the observation data.  These will be used to identify
  // which observations belong to which profile.
  const std::vector<size_t> & record_numbers = obsdb_.recidx_all_recnums();

  oops::Log::debug() << "Unique record numbers" << std::endl;
  for (size_t iProfile : record_numbers)
    oops::Log::debug() << iProfile << ' ';
    oops::Log::debug() << std::endl;

  // Loop over the unique profiles
  for (size_t iProfile : record_numbers) {
    const std::vector<size_t> & obs_numbers = obsdb_.recidx_vector(iProfile);

    // Find the maximum observation of the profile and its index
    float maxObs;
    int maxObsIdx = -1;
    for (size_t iobs : obs_numbers) {
      if (impactParameterObs[0][iobs] > 0 && impactParameterObs[0][iobs] != missingFloat) {
        if (maxObsIdx < 0) {
          maxObs = obsValue[0][iobs];
          maxObsIdx = iobs;
        } else if (obsValue[0][iobs] > obsValue[0][maxObsIdx]) {
          maxObs = std::max(maxObs, obsValue[0][iobs]);
          maxObsIdx = iobs;
        }
      }
    }

    for (size_t iFilterVar = 0; iFilterVar < filtervars.nvars(); ++iFilterVar) {
      float tossMaxHeight = 0;  // initilize the toss height for step 2
      float maxObsVal = 0;  // initilize the maximum obs value for step 2
      for (size_t iobs : obs_numbers) {
        // only check observations below the assigned impact height
        if ((impactParameterObs[0][iobs] - radiusCurvature[0][iobs])
             >= parameters_.checkHeight.value()) continue;
        std::vector<float> refracProfile;
        std::vector<float> heightProfile;
        for (size_t iLevel = 0; iLevel < nRefLevels; ++iLevel)
          if (refractivity[iLevel][iobs] != missingFloat &&
              modelHeights[iLevel][iobs] != missingFloat) {
              refracProfile.push_back(refractivity[iLevel][iobs]);
              heightProfile.push_back(modelHeights[iLevel][iobs]);
          }
        // calculate impact parameter in model space
        const std::vector<float> & impactParameterModel = calcImpactParameterModel(
                        refracProfile, heightProfile,
                        latitude[0][iobs],
                        geoid[0][iobs],
                        radiusCurvature[0][iobs]);
        const std::vector<float> & refGradient = calcRefractivityGradient(
                        refracProfile, impactParameterModel);

        // check from the middle of the model levels to surface
        int ncheckstart = static_cast<int>(refGradient.size()/2.0);
        float gradientThreshold = parameters_.thresholdS1.value()
                                * Constants::superRefractionCritVal;

        for (size_t iLevel = ncheckstart; iLevel < refGradient.size(); ++iLevel) {
          // super refraction check step1
          if (refGradient[iLevel] == missingFloat) {
            flagged[iFilterVar][iobs] = true;
            break;
          }
          float modelImpact = impactParameterModel[iLevel-parameters_.closeLayers.value()];
          if (std::abs(refGradient[iLevel]) >= gradientThreshold &&
            impactParameterObs[0][iobs] <= modelImpact) {
            flagged[iFilterVar][iobs] = true;
            break;
          }
        }  // end the iLevel loop for step1

        // identify the maximum height of super refraction layer for step2 (or 3)
        if (parameters_.step.value() > 1 && obsValue[0][iobs] >= parameters_.maxObs.value()) {
          for (size_t iLevel = ncheckstart; iLevel < refGradient.size(); ++iLevel) {
            if (refGradient[iLevel] >= parameters_.thresholdS2.value() *
                Constants::superRefractionCritVal ) {
               if (obsValue[0][iobs] > maxObsVal && parameters_.step.value() == 2) {
                 maxObsVal = obsValue[0][iobs];
                 tossMaxHeight = impactParameterObs[0][iobs];
               } else {
                 tossMaxHeight = std::max(tossMaxHeight, impactParameterObs[0][iobs]);
               }
               break;
            }
          }  // end the iLevel loop for step2
        }
      }  // end iobs loop

      // apply super refraction check step2 (or 3)
      if (parameters_.step.value() > 1 && tossMaxHeight > 0) {
        for (size_t iobs : obs_numbers) {
          if (impactParameterObs[0][iobs] - radiusCurvature[0][iobs]
              >= parameters_.checkHeight.value()) continue;
          if (impactParameterObs[0][iobs] <= tossMaxHeight) {
            flagged[iFilterVar][iobs] = true;
          }
        }
      }  // end step 2(3) check
    }  // end iFilterVar loop
  }   // end iProfile loop
}  // end applyFilter

std::vector<float> SuperRefractionCheckNBAM::calcImpactParameterModel(
            const std::vector<float> & refrac,
            const std::vector<float> & geopotentialHeight,
            float lat,
            float geoid,
            float radiusCurv) const {
  const float missingFloat = util::missingValue<float>();
  std::vector<float> impactParameterModel;
  float geometricHeight;
  for (size_t iLevel = 0; iLevel < geopotentialHeight.size(); ++iLevel) {
    geometricHeight = formulas::Geopotential_to_Geometric_Height(lat,
                                geopotentialHeight[iLevel]+geoid);
    impactParameterModel.push_back((1.0E-6f*refrac[iLevel]+1)*(radiusCurv+geometricHeight));
  }
  return impactParameterModel;
}

std::vector<float> SuperRefractionCheckNBAM::calcRefractivityGradient(
            const std::vector<float> & refrac,
            const std::vector<float> & height) const {
  const float missingFloat = util::missingValue<float>();
  std::vector<float> gradient;
  for (size_t iLevel = 0; iLevel < height.size(); ++iLevel) {
    if (iLevel > 0 && std::abs(height[iLevel-1] - height[iLevel]) > 0.0) {
      gradient.push_back((refrac[iLevel-1] - refrac[iLevel]) /
                         (height[iLevel-1] - height[iLevel]));
    } else {
      gradient.push_back(missingFloat);
    }
  }
  return gradient;
}

// -----------------------------------------------------------------------------

void SuperRefractionCheckNBAM::print(std::ostream & os) const {
  os << "SuperRefractionCheckNBAM: config = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
