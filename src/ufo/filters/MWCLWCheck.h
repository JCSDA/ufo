/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_MWCLWCHECK_H_
#define UFO_FILTERS_MWCLWCHECK_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/base/ParameterTraitsVariables.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/NumericConstraints.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

/// Parameters controlling the operation of the MWCLWCheck filter.
class MWCLWCheckParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(MWCLWCheckParameters, FilterParametersBase)

 public:
  oops::RequiredParameter<oops::ObsVariables> clwVariables{"clw variables", this};

  oops::RequiredParameter<std::vector<float>> clwThresholds{"clw_thresholds", this};

  /// Controls how the clw is calculated:
  ///
  /// 1. Use observed BTs.
  /// 2. Use calculated BTs.
  /// 3. Symmetric calculation.
  oops::RequiredParameter<int> clwOption{"clw_option", this, {oops::minConstraint<int>(1),
                                                              oops::maxConstraint<int>(3)}};
};

class MWCLWCheck : public FilterBase,
                   private util::ObjectCounter<MWCLWCheck> {
 public:
  /// The type of parameters accepted by the constructor of this filter.
  /// This typedef is used by the FilterFactory.
  typedef MWCLWCheckParameters Parameters_;

  static const std::string classname() {return "ufo::MWCLWCheck";}

  MWCLWCheck(ioda::ObsSpace &, const Parameters_ &,
             std::shared_ptr<ioda::ObsDataVector<int> >,
             std::shared_ptr<ioda::ObsDataVector<float> >);
  ~MWCLWCheck();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::clw;}

  Parameters_ parameters_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_MWCLWCHECK_H_
