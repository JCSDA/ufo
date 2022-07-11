/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_DIFFERENCECHECK_H_
#define UFO_FILTERS_DIFFERENCECHECK_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

/// Parameters controlling the operation of the DifferenceCheck filter.
class DifferenceCheckParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(DifferenceCheckParameters, FilterParametersBase)

 public:
  /// Name of the reference variable.
  oops::RequiredParameter<Variable> ref{"reference", this};
  /// Name of the test variable.
  oops::RequiredParameter<Variable> val{"value", this};

  /// The filter will flag observations for which the difference `test - reference` is below
  /// `minvalue`.
  oops::OptionalParameter<float> minvalue{"minvalue", this};
  /// The filter will flag observations for which the difference `test - reference` is above
  /// `maxvalue`.
  oops::OptionalParameter<float> maxvalue{"maxvalue", this};

  /// If the `threshold` option is specified, the filter behaves as if `minvalue` was set to
  /// `-threshold` and `maxvalue` was set to `threshold` (overriding any values of these options
  /// specified independently).
  oops::OptionalParameter<float> threshold{"threshold", this};
};

/// A filter that compares the difference between a test variable and a reference variable and
/// flags observations for which this difference is outside of a prescribed range.
///
/// See DifferenceCheckParameters for the documentation of the parameters controlling this filter.
class DifferenceCheck : public FilterBase,
                        private util::ObjectCounter<DifferenceCheck> {
 public:
  /// The type of parameters accepted by the constructor of this filter.
  /// This typedef is used by the FilterFactory.
  typedef DifferenceCheckParameters Parameters_;

  static const std::string classname() {return "ufo::DifferenceCheck";}

  DifferenceCheck(ioda::ObsSpace &, const Parameters_ &,
                  std::shared_ptr<ioda::ObsDataVector<int> >,
                  std::shared_ptr<ioda::ObsDataVector<float> >);
  ~DifferenceCheck();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::diffref;}

  Parameters_ parameters_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_DIFFERENCECHECK_H_
