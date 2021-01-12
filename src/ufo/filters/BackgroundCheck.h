/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_BACKGROUNDCHECK_H_
#define UFO_FILTERS_BACKGROUNDCHECK_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/Printable.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

/// Parameters controlling the operation of the BackgroundCheck filter.
class BackgroundCheckParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(BackgroundCheckParameters, FilterParametersBase)

 public:
  /// The filter will flag observations whose bias-corrected value differs from its model equivalent
  /// by more than `threshold` times the current estimate of the observation error.
  ///
  /// `threshold` can be a real number or the name of a variable.
  oops::OptionalParameter<std::string> threshold{"threshold", this};

  /// The filter will flag observations whose bias-corrected value differs from its model equivalent
  /// by more than `absolute threshold`
  ///
  /// `absolute threshold` can be a real number or the name of a variable.
  oops::OptionalParameter<std::string> absoluteThreshold{"absolute threshold", this};

  /// The filter will flag observations whose bias-corrected value differs from its model equivalent
  /// by more than `function absolute threshold`
  ///
  /// `function absolute threshold` should be a list of ObsFunctions (but only the first element
  /// of this list is currently taken into account). The ObsFunctions may return multiple values
  /// per observation, which is especially useful for quality checking data with multiple channels.
  ///
  /// If `function absolute threshold` is set, neither `threshold` nor `absolute threshold` should
  /// be set.
  oops::OptionalParameter<std::vector<Variable>> functionAbsoluteThreshold{
    "function absolute threshold", this};

  /// Name of the HofX group used to replace the default group (default is HofX)
  oops::Parameter<std::string> test_hofx{"test_hofx", "HofX", this};
};

/// BackgroundCheck: check observation closeness to background.
///
/// See BackgroundCheckParameters for the documentation of the parameters controlling this filter.
class BackgroundCheck : public FilterBase,
                        private util::ObjectCounter<BackgroundCheck> {
 public:
  /// The type of parameters accepted by the constructor of this filter.
  /// This typedef is used by the FilterFactory.
  typedef BackgroundCheckParameters Parameters_;

  static const std::string classname() {return "ufo::BackgroundCheck";}

  BackgroundCheck(ioda::ObsSpace &, const Parameters_ &,
                  std::shared_ptr<ioda::ObsDataVector<int> >,
                  std::shared_ptr<ioda::ObsDataVector<float> >);
  ~BackgroundCheck();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::fguess;}

  Parameters_ parameters_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_BACKGROUNDCHECK_H_
