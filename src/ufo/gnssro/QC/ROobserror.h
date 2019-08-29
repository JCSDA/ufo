/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_GNSSRO_QC_ROOBSERROR_H_
#define UFO_GNSSRO_QC_ROOBSERROR_H_

#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "ROobserror.interface.h"

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
  class ObsDiagnostics;

/// ROobserror: calculate observational errors

class ROobserror : public util::Printable,
                   private util::ObjectCounter<ROobserror> {
 public:
  static const std::string classname() {return "ufo::ROobserror";}

  ROobserror(const ioda::ObsSpace &, const eckit::Configuration &,
             boost::shared_ptr<ioda::ObsDataVector<int> >,
             boost::shared_ptr<ioda::ObsDataVector<float> >);
  ~ROobserror();

  void preProcess() const {}
  void priorFilter(const GeoVaLs &) const;
  void postFilter(const ioda::ObsVector &, const ObsDiagnostics &) const {}

  const oops::Variables & requiredGeoVaLs() const {return geovars_;}
  const oops::Variables & requiredHdiagnostics() const {return diagvars_;}

 private:
  void print(std::ostream &) const;
  F90roerr key_;
  const oops::Variables geovars_;
  const oops::Variables diagvars_;
  boost::shared_ptr<ioda::ObsDataVector<int> > flags_;
  boost::shared_ptr<ioda::ObsDataVector<float> > obserr_;
};

}  // namespace ufo

#endif  // UFO_GNSSRO_QC_ROOBSERROR_H_
