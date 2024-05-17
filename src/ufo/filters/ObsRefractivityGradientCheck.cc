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
  const oops::ObsVariables observed = obsdb_.obsvariables();
  const float missingFloat = util::missingValue<float>();

  // Read observation height and atmosphericRefractivity
  Variable obsHeight = Variable("MetaData/height");
  std::vector<float> height;
  data_.get(obsHeight, height);
  Variable obsRefractivity = Variable("ObsValue/atmosphericRefractivity");
  std::vector<float> refractivity;
  data_.get(obsRefractivity, refractivity);

  const std::vector<size_t> & record_numbers = obsdb_.recidx_all_recnums();
  oops::Log::debug() <<"Unique record numbers" << std::endl;
  for (size_t iProfile : record_numbers)
    oops::Log::debug() << iProfile << ' ';
  oops::Log::debug() << std::endl;

  // Loop over the unique profiles
  for (size_t iProfile : record_numbers) {
    const std::vector<size_t> & obs_numbers = obsdb_.recidx_vector(iProfile);

    // Find the set of indices which allow us to sort the variables by height
    // decending
    std::vector<size_t> idx(obs_numbers.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(),
         [&height, &obs_numbers](size_t i1, size_t i2)
         {return height[obs_numbers[i1]] > height[obs_numbers[i2]];});

    std::vector<float> refracProfile;
    std::vector<float> heightProfile;
    for (size_t isort : idx) {
      refracProfile.push_back(refractivity[obs_numbers[isort]]);
      heightProfile.push_back(height[obs_numbers[isort]]);
    }
    const std::vector<float> & gradient = calcVerticalGradient(refracProfile, heightProfile);
    const std::vector<float> & secondDeriv = calcVerticalGradient(gradient, heightProfile);

    for (size_t jv = 0; jv < filtervars.nvars(); ++jv) {
      const size_t iv = observed.find(filtervars.variable(jv).variable());
      for (size_t isort : idx) {
        if ( apply[obs_numbers[isort]] && (*flags_)[iv][obs_numbers[isort]] == QCflags::pass ) {
          if ((gradient[isort] <= parameters_.gradientMin.value() ||
              gradient[isort] >= parameters_.gradientMax.value() ||
              gradient[isort] == 0 ||
              abs(secondDeriv[isort]) >= parameters_.secondDerivative.value()) &&
              (gradient[isort] != missingFloat  && secondDeriv[isort]!= missingFloat) &&
              height[obs_numbers[isort]] < parameters_.maxCheckHeight.value()) {
            // reject all observations below isort
            for (size_t jobs = isort; jobs < idx.size(); ++jobs) {
              if (apply[obs_numbers[jobs]] && (*flags_)[iv][obs_numbers[jobs]] == QCflags::pass) {
                flagged[jv][obs_numbers[jobs]] = true;
              }
            }
            break;
          }
        }
      }  //  end isort loop
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
