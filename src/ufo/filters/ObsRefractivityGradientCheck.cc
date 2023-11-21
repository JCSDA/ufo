/*
 * (C) Copyright 2017-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/ObsRefractivityGradientCheck.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/util/Logger.h"
#include "ufo/filters/QCflags.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsRefractivityGradientCheck::ObsRefractivityGradientCheck(
                                 ioda::ObsSpace & obsdb,
                                 const Parameters_ & parameters,
                                 std::shared_ptr<ioda::ObsDataVector<int> > flags,
                                 std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, parameters, flags, obserr), parameters_(parameters)
{
  oops::Log::trace() << "ObsRefractivityGradientCheck: "
                     << "using GNSSRO super refraction check method" << std::endl;
  allvars_ += Variable("ObsValue/atmosphericRefractivity");
  allvars_ += Variable("MetaData/height");
  if (obsdb.obs_group_vars().empty())
    throw eckit::BadParameter("group variables configuration is empty.", Here());
}
// -----------------------------------------------------------------------------

ObsRefractivityGradientCheck::~ObsRefractivityGradientCheck() {
  oops::Log::trace() << "ObsRefractivityGradientCheck: destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRefractivityGradientCheck::applyFilter(
                                      const std::vector<bool> & apply,
                                      const Variables & filtervars,
                                      std::vector<std::vector<bool>> & flagged) const {
  oops::Log::trace() << "ObsRefractivityGradientCheck postFilter" << std::endl;
  const oops::Variables observed = obsdb_.obsvariables();
  const float missingFloat = util::missingValue<float>();

  ioda::ObsDataVector<float> refractivity(obsdb_, "atmosphericRefractivity", "ObsValue");
  ioda::ObsDataVector<float> height(obsdb_, "height", "MetaData");

  const std::vector<size_t> & record_numbers = obsdb_.recidx_all_recnums();
  oops::Log::debug() <<"Unique record numbers" << std::endl;
  for (size_t iProfile : record_numbers)
    oops::Log::debug() << iProfile << ' ';
  oops::Log::debug() << std::endl;

  // Loop over the unique profiles
  for (size_t iProfile : record_numbers) {
    const std::vector<size_t> & obs_numbers = obsdb_.recidx_vector(iProfile);
    std::vector<float> refracProfile;
    std::vector<float> heightProfile;

    for (size_t iobs : obs_numbers) {
      refracProfile.push_back(refractivity[0][iobs]);
      heightProfile.push_back(height[0][iobs]);
    }

    const std::vector<float> & gradient = calcVerticalGradient(refracProfile, heightProfile);
    const std::vector<float> & secondDeriv = calcVerticalGradient(gradient, heightProfile);

    refracProfile.clear();
    heightProfile.clear();

    for (size_t jv = 0; jv < filtervars.nvars(); ++jv) {
      int jobs = 0;
      const size_t iv = observed.find(filtervars.variable(jv).variable());
      for (size_t iobs : obs_numbers) {
        if ( apply[iobs] && (*flags_)[iv][iobs] == QCflags::pass ) {
          if ((gradient[jobs] <= parameters_.gradientMin.value() ||
              gradient[jobs] >= parameters_.gradientMax.value() ||
              gradient[jobs] == 0 ||
              abs(secondDeriv[jobs]) >= parameters_.secondDerivative.value()) &&
              (gradient[jobs] != missingFloat  && secondDeriv[jobs]!= missingFloat)) {
            flagged[jv][iobs] = true;
          }
        }
       jobs++;
      }
    }  //  end jv loop
  }  //  end iProfile loop
}

std::vector<float> ObsRefractivityGradientCheck::calcVerticalGradient(
        const std::vector<float> & refrac,
        const std::vector<float> & height) const {
  std::vector<float> gradient;
  const float missingFloat = util::missingValue<float>();
  const float hugeValue = 1.0e8;

  gradient.push_back(missingFloat);
  for (size_t iLevel = 1; iLevel < height.size(); ++iLevel) {
    if ( height[iLevel] != missingFloat && refrac[iLevel] != missingFloat &&
         height[iLevel-1] != missingFloat && refrac[iLevel-1] != missingFloat &&
         refrac[iLevel] < hugeValue && refrac[iLevel-1] < hugeValue) {
      gradient.push_back((refrac[iLevel-1] - refrac[iLevel]) / (height[iLevel-1] - height[iLevel]));
    } else {
      gradient.push_back(missingFloat);
    }
  }

  return gradient;
}

// -----------------------------------------------------------------------------

void ObsRefractivityGradientCheck::print(std::ostream & os) const {
  os << "ObsRefractivityGradientCheck: config = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
