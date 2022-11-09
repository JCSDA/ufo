/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_GNSSRO_QC_ROOBSERROR_H_
#define UFO_OPERATORS_GNSSRO_QC_ROOBSERROR_H_

#include <Eigen/Dense>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/filters/FilterBase.h"
#include "ROobserror.interface.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

/// ROobserror: calculate observational errors

class ROobserror : public FilterBase,
                   private util::ObjectCounter<ROobserror> {
 public:
  static const std::string classname() {return "ufo::ROobserror";}

  ROobserror(ioda::ObsSpace &, const eckit::Configuration &,
             std::shared_ptr<ioda::ObsDataVector<int> >,
             std::shared_ptr<ioda::ObsDataVector<float> >);
  ~ROobserror();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return 76;}
  int n_horiz = 1;

  F90roerr key_;
};

}  // namespace ufo

#endif  // UFO_OPERATORS_GNSSRO_QC_ROOBSERROR_H_
