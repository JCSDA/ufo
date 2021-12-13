/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_OBSBOUNDSCHECK_H_
#define UFO_FILTERS_OBSBOUNDSCHECK_H_

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

/// Parameters controlling the operation of the ObsBoundsCheck filter.
class ObsBoundsCheckParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsBoundsCheckParameters, FilterParametersBase)

 public:
  /// Minimum allowed value of the tested variables.
  oops::OptionalParameter<float> minvalue{"minvalue", this};

  /// Maximum allowed value of the tested variables.
  oops::OptionalParameter<float> maxvalue{"maxvalue", this};

  /// Variables to be compared against the bounds specified in the `minvalue` and `maxvalue`
  /// options.
  ///
  /// There are three valid possibilities:
  ///
  /// * If this option is not set, the filter will flag each filter variable at each location
  ///   where the measured value of that variable lies outside the specified bounds.
  ///
  /// * If this option is set to a single-element list containing only one single-channel variable
  ///   or the `flag all filter variables if any test variable is out of bounds` option is set to
  ///   `true`, the filter will flag each filter variable at each location where any test variable
  ///   lies outside the specified bounds. If
  ///   `test only filter variables with passed qc when flagging all filter variables` is set to
  ///   true and the number of filter vars is not equal to the number of test vars an error
  ///   will be thrown.
  ///
  /// * If this option is set to a list with as many elements as there are filter variables and
  ///   the `flag all filter variables if any test variable is out of bounds` option is set to
  ///   `false`, the filter will flag each filter variable at each location
  ///   where the corresponding test variable lies outside the specified bounds.
  oops::OptionalParameter<std::vector<Variable>> testVariables{"test variables", this};

  /// Set this option to `true` to flag all filter variables at each location where any test
  /// variable lies outside the specified bounds.
  ///
  /// This option is ignored if the `test variables` option is not set.
  oops::Parameter<bool> flagAllFilterVarsIfAnyTestVarOutOfBounds{
    "flag all filter variables if any test variable is out of bounds", false, this};

  /// Set this option to `true` to only test the current filter variable if it is flagged
  /// as pass qc when using "flag all filter variables if any test variables out of bounds".
  ///
  /// This option is not used if "flagAllFilterVarsIfAnyTestVarOutOfBounds" is false.
  oops::Parameter<bool> onlyTestGoodFilterVarsForFlagAllFilterVars{
    "test only filter variables with passed qc when flagging all filter variables", false, this};

  /// By default, the filter flags filter variables at locations where the corresponding test
  /// variable is set to the missing value indicator. Set this option to `false` to stop it from
  /// doing so (hence assuming "optimistically" that the test variable was in fact in bounds).
  oops::Parameter<bool> treatMissingAsOutOfBounds{"treat missing as out of bounds", true, this};
};

/// \brief Flag observations that lie outside specified bounds.
///
/// This is a generic quality control filter based on observation data only.
///
/// See ObsBoundsCheckParameters for the documentation of the parameters controlling this filter.

class ObsBoundsCheck : public FilterBase,
                       private util::ObjectCounter<ObsBoundsCheck> {
 public:
  /// The type of parameters accepted by the constructor of this filter.
  /// This typedef is used by the FilterFactory.
  typedef ObsBoundsCheckParameters Parameters_;

  static const std::string classname() {return "ufo::ObsBoundsCheck";}

  ObsBoundsCheck(ioda::ObsSpace &, const Parameters_ &,
                 std::shared_ptr<ioda::ObsDataVector<int> >,
                 std::shared_ptr<ioda::ObsDataVector<float> >);
  ~ObsBoundsCheck();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::bounds;}
  Parameters_ parameters_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_OBSBOUNDSCHECK_H_
