/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_PREQC_H_
#define UFO_FILTERS_PREQC_H_

#include <ostream>

#include "boost/shared_ptr.hpp"

#include "eckit/config/LocalConfiguration.h"
#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/Printable.h"

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
class GeoVaLs;

class PreQC : public util::Printable {
 public:
  PreQC(ioda::ObsSpace &, const eckit::Configuration &,
        boost::shared_ptr<ioda::ObsDataVector<int> >,
        boost::shared_ptr<ioda::ObsDataVector<float> >);
  ~PreQC() {}

  void priorFilter(const GeoVaLs &) const {}
  void postFilter(const ioda::ObsVector &) const {}

  const oops::Variables & requiredGeoVaLs() const {return nogeovals_;}

 private:
  void print(std::ostream &) const;

  const oops::Variables nogeovals_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_PREQC_H_
