/* * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/BackgroundCheck.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <string>
#include <type_traits>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/util/Logger.h"

#include "ufo/filters/getScalarOrFilterData.h"
#include "ufo/filters/QCflags.h"

namespace ufo {

// -----------------------------------------------------------------------------

BackgroundCheck::BackgroundCheck(ioda::ObsSpace & obsdb, const Parameters_ & parameters,
                                 std::shared_ptr<ioda::ObsDataVector<int> > flags,
                                 std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, parameters, flags, obserr, VariableNameMap(parameters.AliasFile.value())),
    parameters_(parameters)
{
  oops::Log::trace() << "BackgroundCheck constructor" << std::endl;

  // Typical use would be HofX group, but during testing, we include option for GsiHofX
  std::string test_hofx = parameters_.test_hofx.value();

  if (parameters_.functionAbsoluteThreshold.value()) {
    for (const Variable & var : *(parameters_.functionAbsoluteThreshold.value()))
      allvars_ += var;
  }
  allvars_ += Variables(filtervars_, test_hofx);
  // Add BG error variable if threshold is required wrt BG error:
  for (size_t jv = 0; jv < filtervars_.size(); ++jv) {
    if (parameters_.thresholdWrtBGerror.value()) {
      allvars_ += backgrErrVariable(filtervars_[jv]);
    }
  }
  ASSERT(parameters_.threshold.value() ||
         parameters_.absoluteThreshold.value() ||
         parameters_.functionAbsoluteThreshold.value());
  if (parameters_.functionAbsoluteThreshold.value()) {
    ASSERT(!parameters_.threshold.value() &&
           !parameters_.absoluteThreshold.value());
    ASSERT(!parameters_.functionAbsoluteThreshold.value()->empty());
  }
  if (parameters_.thresholdWrtBGerror.value()) {
    ASSERT(parameters_.threshold.value());
  }
}

// -----------------------------------------------------------------------------

BackgroundCheck::~BackgroundCheck() {
  oops::Log::trace() << "BackgroundCheck destructed" << std::endl;
}

// -----------------------------------------------------------------------------
/// Return the name of the variable containing the background error estimate of the
/// specified filter variable.

Variable BackgroundCheck::backgrErrVariable(const Variable &filterVariable) const {
  return Variable("ObsDiag/" +
                  nameMap_.convertName(filterVariable.variable()).name() +
                  "_background_error");
}

// -----------------------------------------------------------------------------

void BackgroundCheck::applyFilter(const std::vector<bool> & apply,
                                  const Variables & filtervars,
                                  std::vector<std::vector<bool>> & flagged) const {
  oops::Log::trace() << "BackgroundCheck postFilter" << std::endl;
  const oops::ObsVariables observed = obsdb_.obsvariables();
  const float missing = util::missingValue<float>();
  oops::Log::debug() << "BackgroundCheck obserr: " << *obserr_ << std::endl;
  ioda::ObsDataVector<float> obs(obsdb_, filtervars.toOopsObsVariables(), "ObsValue");

  std::string test_hofx = parameters_.test_hofx.value();
  Variables varhofx(filtervars, test_hofx);

// Get function absolute threshold
  if (parameters_.functionAbsoluteThreshold.value()) {
//  Get function absolute threshold info from configuration
    const Variable &rtvar = parameters_.functionAbsoluteThreshold.value()->front();
    ioda::ObsDataVector<float> function_abs_threshold(obsdb_, rtvar.toOopsObsVariables());
    data_.get(rtvar, function_abs_threshold);

    for (size_t jv = 0; jv < filtervars.nvars(); ++jv) {
      size_t iv = observed.find(filtervars.variable(jv).variable());
      //    H(x)
      std::vector<float> hofx;
      data_.get(varhofx.variable(jv), hofx);
      for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
        // style note: use missingValue<decltype(foo)> because obserr_ is declared in another file
        // and its underlying type could in principle change without it being obvious here.
        const bool obserrNotMissing = ((*obserr_)[iv][jobs]
                                       != util::missingValue<std::remove_reference_t<
                                                             decltype((*obserr_)[iv][jobs])>>());
        if (apply[jobs] && (*flags_)[iv][jobs] == QCflags::pass && obserrNotMissing) {
          ASSERT(obs[jv][jobs] != missing);
          ASSERT(hofx[jobs] != missing);

//        Threshold for current observation
          float zz = function_abs_threshold[jv][jobs];
//        Check distance from background
          if (std::abs(static_cast<float>(hofx[jobs]) - obs[jv][jobs]) > zz) {
            flagged[jv][jobs] = true;
          }
        }
      }
    }
  } else {
    Variables varbias(filtervars, "ObsBiasData");
    for (size_t jv = 0; jv < filtervars.nvars(); ++jv) {
      size_t iv = observed.find(filtervars.variable(jv).variable());
//    H(x) (including bias correction)
      std::vector<float> hofx;
      data_.get(varhofx.variable(jv), hofx);
//    H(x) error
      std::vector<float> hofxerr;
      bool thresholdWrtBGerror = parameters_.thresholdWrtBGerror.value();
      if (thresholdWrtBGerror) {
        data_.get(backgrErrVariable(filtervars[jv]), hofxerr);
      }
//    Bias correction (only read in if removeBiasCorrection is set to true, otherwise
//    set to zero).
      std::vector<float> bias(obsdb_.nlocs(), 0.0);
      if (parameters_.removeBiasCorrection) {
        data_.get(varbias.variable(jv), bias);
      }

//    Threshold for current variable
      std::vector<float> abs_thr(obsdb_.nlocs(), std::numeric_limits<float>::max());
      std::vector<float> thr(obsdb_.nlocs(), std::numeric_limits<float>::max());
      if (parameters_.absoluteThreshold.value())
        abs_thr = getScalarOrFilterData(*parameters_.absoluteThreshold.value(), data_);
      if (parameters_.threshold.value())
        thr = getScalarOrFilterData(*parameters_.threshold.value(), data_);

      for (size_t jobs = 0; jobs < obsdb_.nlocs(); ++jobs) {
        // style note: use missingValue<decltype(foo)> because obserr_ is declared in another file
        // and its underlying type could in principle change without it being obvious here.
        const bool obserrNotMissing = ((*obserr_)[iv][jobs]
                                       != util::missingValue<std::remove_reference_t<
                                                             decltype((*obserr_)[iv][jobs])>>());
        if (apply[jobs] && (*flags_)[iv][jobs] == QCflags::pass && obserrNotMissing) {
          ASSERT(obs[jv][jobs] != missing);
          ASSERT(hofx[jobs] != missing);
          ASSERT(bias[jobs] != missing);

          const std::vector<float> &errorMultiplier = thresholdWrtBGerror ?
                                                      hofxerr : (*obserr_)[iv];
//        Threshold for current observation
          float zz = (thr[jobs] == std::numeric_limits<float>::max()) ? abs_thr[jobs] :
            std::min(abs_thr[jobs], thr[jobs] * errorMultiplier[jobs]);

          ASSERT(zz < std::numeric_limits<float>::max() && zz > 0.0);

//        Check distance from background. hofx includes bias correction.
//        If removeBiasCorrection is set to true, `bias` contains bias correction, and
//           it is removed from hofx.
//        Otherwise, `bias` is set to zero, and bias correction is not removed from hofx.
          if (std::abs(hofx[jobs] - obs[jv][jobs] - bias[jobs]) > zz) {
            flagged[jv][jobs] = true;
          }
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------

void BackgroundCheck::print(std::ostream & os) const {
  os << "BackgroundCheck::print not yet implemented ";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
