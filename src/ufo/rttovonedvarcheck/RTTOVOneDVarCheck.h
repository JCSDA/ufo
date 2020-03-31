/*
 * (C) Copyright 2017-2020 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_RTTOVONEDVARCHECK_RTTOVONEDVARCHECK_H_
#define UFO_RTTOVONEDVARCHECK_RTTOVONEDVARCHECK_H_

#include <cmath>
#include <ostream>
#include <string>
#include <vector>

#include "boost/shared_ptr.hpp"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"
#include "ufo/rttovonedvarcheck/RTTOVOneDVarCheck.interface.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

/// RTTOVOneDVarCheck: check observation closeness to background

class RTTOVOneDVarCheck : public FilterBase,
                     private util::ObjectCounter<RTTOVOneDVarCheck> {
 public:
  static const std::string classname() {return "ufo::RTTOVOneDVarCheck";}

  RTTOVOneDVarCheck(ioda::ObsSpace &, const eckit::Configuration &,
                  boost::shared_ptr<ioda::ObsDataVector<int> >,
                  boost::shared_ptr<ioda::ObsDataVector<float> >);
  ~RTTOVOneDVarCheck();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::onedvar;}

  F90onedvarcheck key_;
  const eckit::LocalConfiguration config_;
  std::vector<int> channels_;
};

}  // namespace ufo

#endif  // UFO_RTTOVONEDVARCHECK_RTTOVONEDVARCHECK_H_
