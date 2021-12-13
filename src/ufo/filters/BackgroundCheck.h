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
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

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
  /// by more than `threshold` times the current estimate of the observation error. Or if the
  /// option "threshold wrt background error" is true, the `threshold` is multiplied by the
  /// background error rather than observation error. E.g.
  ///
  ///   filter variables:
  ///   - name: sea_surface_height
  ///   threshold wrt background error: true
  ///   threshold: 3.0
  ///
  /// `threshold` can be a real number or the name of a variable.
  oops::OptionalParameter<std::string> threshold{"threshold", this};

  /// A switch indicating whether threshold must be multiplied by background error rather than
  /// observation error. If true, `threshold` must have a value.
  oops::Parameter<bool> thresholdWrtBGerror{"threshold wrt background error", false, this};

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

  /// The filter uses bias-corrected H(x) unless `remove bias correction` is set to true.
  oops::Parameter<bool> removeBiasCorrection{"remove bias correction", false, this};

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
  /// \brief Return the name of the variable containing the background error estimate of the
  /// specified filter variable.
  Variable backgrErrVariable(const Variable & filterVariable) const;

  Parameters_ parameters_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_BACKGROUNDCHECK_H_
