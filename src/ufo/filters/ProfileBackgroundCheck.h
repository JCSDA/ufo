/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_PROFILEBACKGROUNDCHECK_H_
#define UFO_FILTERS_PROFILEBACKGROUNDCHECK_H_

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

/// Parameters controlling the operation of the ProfileBackgroundCheck filter.
class ProfileBackgroundCheckParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(ProfileBackgroundCheckParameters, FilterParametersBase)

 public:
  /// The filter will flag observations whose bias-corrected value differs from its model equivalent
  /// by more than `relative threshold` times the current estimate of the observation error.
  ///
  /// `relative threshold` can be a real number or the name of a variable.
  oops::OptionalParameter<std::string> relativeThreshold{"relative threshold", this};

  /// The filter will flag observations whose bias-corrected value differs from its model equivalent
  /// by more than `absolute threshold`
  ///
  /// `absolute threshold` can be a real number or the name of a variable.
  oops::OptionalParameter<std::string> absoluteThreshold{"absolute threshold", this};
};

/// ProfileBackgroundCheck: check observation closeness to background over a profile
///
/// See the cc file for more details.

class ProfileBackgroundCheck : public FilterBase,
                        private util::ObjectCounter<ProfileBackgroundCheck> {
 public:
  /// The type of parameters accepted by the constructor of this filter.
  /// This typedef is used by the FilterFactory.
  typedef ProfileBackgroundCheckParameters Parameters_;

  static const std::string classname() {return "ufo::ProfileBackgroundCheck";}

  ProfileBackgroundCheck(ioda::ObsSpace &, const Parameters_ &,
                  std::shared_ptr<ioda::ObsDataVector<int> >,
                  std::shared_ptr<ioda::ObsDataVector<float> >);
  ~ProfileBackgroundCheck();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::fguess;}
  Parameters_ parameters_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_PROFILEBACKGROUNDCHECK_H_
