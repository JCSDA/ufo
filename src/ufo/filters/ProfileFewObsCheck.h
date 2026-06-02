/*
 * (C) British Crown Copyright 2021 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_PROFILEFEWOBSCHECK_H_
#define UFO_FILTERS_PROFILEFEWOBSCHECK_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/NumericConstraints.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

/// Parameters controlling the operation of the ProfileFewObsCheck filter.
class ProfileFewObsCheckParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(ProfileFewObsCheckParameters, FilterParametersBase)

 public:
  /// The filter will flag profiles which contain fewer than `threshold` number
  /// of observations
  oops::OptionalParameter<int> threshold{"threshold", this, {oops::minConstraint(0)}};
  oops::OptionalParameter<float> fraction{"fraction", this,
    {oops::minConstraint<float>(0), oops::maxConstraint<float>(1)}};
};

/// ProfileFewObsCheck: Check the number of observations in a profile
///
/// See the cc file for more details.

class ProfileFewObsCheck : public FilterBase,
                           private util::ObjectCounter<ProfileFewObsCheck> {
 public:
  /// The type of parameters accepted by the constructor of this filter.
  /// This typedef is used by the FilterFactory.
  typedef ProfileFewObsCheckParameters Parameters_;

  static const std::string classname() {return "ufo::ProfileFewObsCheck";}

  ProfileFewObsCheck(ioda::ObsSpace &, const Parameters_ &,
                     std::shared_ptr<ioda::ObsDataVector<int> >,
                     std::shared_ptr<ioda::ObsDataVector<float> >);
  ~ProfileFewObsCheck();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::profile;}
  Parameters_ parameters_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_PROFILEFEWOBSCHECK_H_
