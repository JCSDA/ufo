/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_GNSSRO_QC_BACKGROUNDCHECKROGSI_H_
#define UFO_GNSSRO_QC_BACKGROUNDCHECKROGSI_H_

#include <ostream>
#include <string>

#include "boost/shared_ptr.hpp"

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

/// BackgroundCheckROGSI: check observation closeness to background

class BackgroundCheckROGSI : public util::Printable,
                        private util::ObjectCounter<BackgroundCheckROGSI> {
 public:
  static const std::string classname() {return "ufo::BackgroundCheckROGSI";}

  BackgroundCheckROGSI(ioda::ObsSpace &, const eckit::Configuration &,
                       boost::shared_ptr<ioda::ObsDataVector<int> >,
                       boost::shared_ptr<ioda::ObsDataVector<float> >);
  ~BackgroundCheckROGSI();

  void preProcess() const {}
  void priorFilter(const GeoVaLs &) const;
  void postFilter(const ioda::ObsVector &, const ObsDiagnostics &) const;

  const oops::Variables & requiredGeoVaLs() const {return geovars_;}
  const oops::Variables & requiredHdiagnostics() const {return diagvars_;}

 private:
  void print(std::ostream &) const;

  ioda::ObsSpace & obsdb_;
  const eckit::LocalConfiguration config_;
  const GeoVaLs mutable * gv_;
  const oops::Variables geovars_;
  const oops::Variables diagvars_;
  ioda::ObsDataVector<int> & flags_;
};

}  // namespace ufo

#endif  // UFO_GNSSRO_QC_BACKGROUNDCHECKROGSI_H_
