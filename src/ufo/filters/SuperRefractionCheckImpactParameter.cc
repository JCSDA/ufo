/*
 * (C) Copyright 2017-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/SuperRefractionCheckImpactParameter.h"
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

SuperRefractionCheckImpactParameter::SuperRefractionCheckImpactParameter(
                                 ioda::ObsSpace & obsdb,
                                 const Parameters_ & parameters,
                                 std::shared_ptr<ioda::ObsDataVector<int> > flags,
                                 std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, parameters, flags, obserr), parameters_(parameters)
{
  oops::Log::trace() << "SuperRefractionCheckImpactParameter"<< std::endl;
  allvars_ += Variable("ObsDiag/atmosphericRefractivity_model");
  allvars_ += Variable("ObsDiag/geopotentialHeight_model");
  allvars_ += Variable("MetaData/impactParameterRO");
  allvars_ += Variable("MetaData/earthRadiusCurvature");
  allvars_ += Variable("MetaData/latitude");
  allvars_ += Variable("MetaData/geoidUndulation");

  if (obsdb.obs_group_vars().empty())
    throw eckit::BadParameter("group variables configuration is empty.", Here());
}
// -----------------------------------------------------------------------------

SuperRefractionCheckImpactParameter::~SuperRefractionCheckImpactParameter() {
  oops::Log::trace() << "SuperRefractionCheckImpactParameter: destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void SuperRefractionCheckImpactParameter::applyFilter(
                                      const std::vector<bool> & apply,
                                      const Variables & filtervars,
                                      std::vector<std::vector<bool>> & flagged) const {
  oops::Log::trace() << "SuperRefractionCheckImpactParameter postFilter, "
                     << "using profile check = "<< parameters_.profileCheck.value()<< std::endl;
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
  oops::Log::trace() << data_.nlevs(modelHeightsVariable) << std::endl;
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

    //  use geovals of the bottom observation of the profile if profileCheck=true
    if (parameters_.profileCheck.value()) {
      // Find the observation with the smallest impact parameter (i.e. the
      // bottom of the profile).  If all impactParameters are missing, then this
      // should choose the first observation.
      int bottomOb = -1;
      for (size_t iobs : obs_numbers) {
        if (impactParameterObs[0][iobs] > 0 && impactParameterObs[0][iobs] != missingFloat) {
          if (bottomOb < 0) {
            bottomOb = iobs;
          } else if (impactParameterObs[0][iobs] < impactParameterObs[0][bottomOb]) {
            bottomOb = iobs;
          }
        }
      }
      std::vector<float> refracProfile;
      std::vector<float> heightProfile;

      // count the number of valid levels of model profiles for iobs
      int nLevel = 0;
      for (size_t iLevel = 0; iLevel < nRefLevels; ++iLevel)
         if (refractivity[iLevel][bottomOb] != missingFloat &&
             modelHeights[iLevel][bottomOb] != missingFloat) {
             refracProfile.push_back(refractivity[iLevel][bottomOb]);
             heightProfile.push_back(modelHeights[iLevel][bottomOb]);
             nLevel++;
         }
      // skip to next observation if model profiles are missing for iobs
      if (nLevel == 0) continue;

      const  std::vector<float>  impactParameterModel = calcImpactParameterModel(
                        refracProfile, heightProfile,
                        latitude[0][bottomOb],
                        geoid[0][bottomOb],
                        radiusCurvature[0][bottomOb]);
      float impactParaDiff = missingFloat;
      int kLevel = 0;
      for (size_t iLevel = 1; iLevel < nRefLevels; ++iLevel) {
        impactParaDiff = impactParameterModel[iLevel-1] - impactParameterModel[iLevel];
        if (impactParaDiff <= parameters_.threshold.value()) {
          kLevel = iLevel-1;
          break;
        }
      }
      for (size_t iVar = 0; iVar < filtervars.nvars(); ++iVar) {
        for (size_t iobs : obs_numbers) {
          if (apply[iobs] &&  (*flags_)[iVar][iobs] == QCflags::pass &&
              impactParameterObs[0][iobs] <= impactParameterModel[kLevel] &&
              kLevel > 0) {
            flagged[iVar][iobs] = true;
          }
        }  // end iobs loop
      }  // end iVar loop

    } else {  //  profileCheck = false
      for (size_t iVar = 0; iVar < filtervars.nvars(); ++iVar) {
        for (size_t iobs : obs_numbers) {
          if (apply[iobs] &&  (*flags_)[iVar][iobs] == QCflags::pass) {
            std::vector<float> refracProfile;
            std::vector<float> heightProfile;
            // count the number of valid levels of model profiles for iobs
            int nLevel = 0;
            for (size_t iLevel = 0; iLevel < nRefLevels; ++iLevel)
              if (refractivity[iLevel][iobs] != missingFloat &&
                  modelHeights[iLevel][iobs] != missingFloat) {
                refracProfile.push_back(refractivity[iLevel][iobs]);
                heightProfile.push_back(modelHeights[iLevel][iobs]);
                nLevel++;
              }
            // skip to next observation if model profiles are missing for iobs
            if (nLevel == 0) continue;

            const  std::vector<float>  impactParameterModel = calcImpactParameterModel(
                          refracProfile, heightProfile,
                          latitude[0][iobs],
                          geoid[0][iobs],
                          radiusCurvature[0][iobs]);
            float impactParaDiff = missingFloat;
            int kLevel = 0;
            for (size_t iLevel = 1; iLevel < nRefLevels; ++iLevel) {
              impactParaDiff = impactParameterModel[iLevel-1] - impactParameterModel[iLevel];
              if (impactParaDiff <= parameters_.threshold.value()) {
                kLevel = iLevel-1;
                break;
              }
            }
            if (impactParameterObs[0][iobs] <= impactParameterModel[kLevel] &&
                kLevel > 0)
              flagged[iVar][iobs] = true;
          }
        }  // end iobs loop
      }  // end iVar loop
    }  // end if profileCheck
  }   //  end iProfile loop
}  // end applyFilter

std::vector<float> SuperRefractionCheckImpactParameter::calcImpactParameterModel(
            const std::vector<float> & refrac,
            const std::vector<float> & geopotentialHeight,
            float lat,
            float geoid,
            float radiusCurv) const {
  const float missingFloat = util::missingValue<float>();
  std::vector<float> impactParameterModel;
  float geometricHeight;
  float temp;
  for (size_t iLevel = 0; iLevel < geopotentialHeight.size(); ++iLevel) {
    geometricHeight = formulas::Geopotential_to_Geometric_Height(lat,
                                geopotentialHeight[iLevel]+geoid);
    impactParameterModel.push_back((1.0E-6f*refrac[iLevel]+1)*(radiusCurv+geometricHeight));
  }
  return impactParameterModel;
}

// -----------------------------------------------------------------------------

void SuperRefractionCheckImpactParameter::print(std::ostream & os) const {
  os << "SuperRefractionCheckImpactParameter: config = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
