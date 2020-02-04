/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_ONEDVARCHECK_H_
#define UFO_FILTERS_ONEDVARCHECK_H_

#include <ostream>
#include <string>
#include <vector>
#include <cmath>

#include "boost/shared_ptr.hpp"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"
#include "ufo/onedvarcheck/OneDVarFortran.interface.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

/// OneDVarCheck: check observation closeness to background

class OneDVarCheck : public FilterBase,
                     private util::ObjectCounter<OneDVarCheck> {
 public:
  static const std::string classname() {return "ufo::OneDVarCheck";}

  OneDVarCheck(ioda::ObsSpace &, const eckit::Configuration &,
                  boost::shared_ptr<ioda::ObsDataVector<int> >,
                  boost::shared_ptr<ioda::ObsDataVector<float> >);
  ~OneDVarCheck();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::onedvar;}

  F90onedvarcheck key_;
  const eckit::LocalConfiguration config_;
  std::vector<int> channels_;
  const GeoVaLs mutable * gv_;
  double cost_converge_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_ONEDVARCHECK_H_
