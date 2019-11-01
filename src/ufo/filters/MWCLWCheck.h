/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_MWCLWCHECK_H_
#define UFO_FILTERS_MWCLWCHECK_H_

#include <ostream>
#include <string>
#include <vector>

#include "boost/shared_ptr.hpp"

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
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

/// MWCLWCheck: generic quality control based on observation data only

// Check that observations are within some bounds over some domain

class MWCLWCheck : public FilterBase,
                   private util::ObjectCounter<MWCLWCheck> {
 public:
  static const std::string classname() {return "ufo::MWCLWCheck";}

  MWCLWCheck(ioda::ObsSpace &, const eckit::Configuration &,
             boost::shared_ptr<ioda::ObsDataVector<int> >,
             boost::shared_ptr<ioda::ObsDataVector<float> >);
  ~MWCLWCheck();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const oops::Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::clw;}

  oops::Variables invars_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_MWCLWCHECK_H_
