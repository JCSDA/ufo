/*
 * (C) Copyright 2018-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSDERIVATIVECHECK_H_
#define UFO_FILTERS_OBSDERIVATIVECHECK_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/missingValues.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/Parameter.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

/// Parameters controlling the operation of the ObsDerivativeCheck filter.
class ObsDerivativeCheckParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsDerivativeCheckParameters, FilterParametersBase)

 public:
  oops::Parameter<std::string> independent
    {"independent",
     "Name of the independent variable",
     "",
     this};

  oops::Parameter<std::string> dependent
    {"dependent",
     "Name of the dependent variable",
     "",
     this};

  oops::Parameter<size_t> i1
    {"i1",
     "First index to use for the derivative computation. "
     "If both i1 and i2 are zero then a local derivative will be computed.",
     0,
     this};

  oops::Parameter<size_t> i2
    {"i2",
     "Second index to use for the derivative computation. "
     "If both i1 and i2 are zero then a local derivative will be computed.",
     0,
     this};

  oops::Parameter<float> minvalue
    {"minvalue",
     "An observation will be flagged if its derivative is lower than this bound.",
     util::missingValue<float>(),
     this};

  oops::Parameter<float> maxvalue
    {"maxvalue",
     "An observation will be flagged if its derivative is larger than this bound.",
     util::missingValue<float>(),
     this};
};

/// Derivative check: check if the derivative of one variable with respect to another
/// is within some range.
///
/// See ObsDerivativeCheckParameters for the documentation of the parameters controlling
/// this filter.
class ObsDerivativeCheck : public FilterBase,
                           private util::ObjectCounter<ObsDerivativeCheck> {
 public:
  typedef ObsDerivativeCheckParameters Parameters_;

  static const std::string classname() {return "ufo::ObsDerivativeCheck";}

  ObsDerivativeCheck(ioda::ObsSpace &, const Parameters_ &,
                     std::shared_ptr<ioda::ObsDataVector<int> >,
                     std::shared_ptr<ioda::ObsDataVector<float> >);
  ~ObsDerivativeCheck();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::derivative;}

  Parameters_ parameters_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_OBSDERIVATIVECHECK_H_
