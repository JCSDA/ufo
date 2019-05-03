/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_GNSSRO_QC_ROGEOREALITYCHECK_H_
#define UFO_GNSSRO_QC_ROGEOREALITYCHECK_H_

#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "ROgeorealityCheck.interface.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;

/// ROgeorealityCheck:RO geophysical reality check

class ROgeorealityCheck : public util::Printable,
                        private util::ObjectCounter<ROgeorealityCheck> {
 public:
  static const std::string classname() {return "ufo::ROgeorealityCheck";}

  ROgeorealityCheck(const ioda::ObsSpace &, const eckit::Configuration &,
                    boost::shared_ptr<ioda::ObsDataVector<int> >,
                    boost::shared_ptr<ioda::ObsDataVector<float> >);
  ~ROgeorealityCheck();

  void priorFilter(const GeoVaLs &) const;
  void postFilter(const ioda::ObsVector &) const;

  const oops::Variables & requiredGeoVaLs() const {return nogeovals_;}

 private:
  void print(std::ostream &) const;
  F90rogeorealitycheck key_;
  const oops::Variables nogeovals_;
  boost::shared_ptr<ioda::ObsDataVector<int> > flags_;
};

}  // namespace ufo

#endif  // UFO_GNSSRO_QC_ROGEOREALITYCHECK_H_
