/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_OBSDIAGNOSTICSWRITER_H_
#define UFO_FILTERS_OBSDIAGNOSTICSWRITER_H_

#include <ostream>
#include <string>

#include "boost/shared_ptr.hpp"

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/Printable.h"
#include "ufo/filters/Variables.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
class GeoVaLs;
class ObsDiagnostics;

class ObsDiagnosticsWriter : public util::Printable,
                             private util::ObjectCounter<ObsDiagnosticsWriter> {
 public:
  static const std::string classname() {return "ufo::ObsDiagnosticsWriter";}

  ObsDiagnosticsWriter(ioda::ObsSpace &, const eckit::Configuration &,
                       boost::shared_ptr<ioda::ObsDataVector<int> >,
                       boost::shared_ptr<ioda::ObsDataVector<float> >);
  ~ObsDiagnosticsWriter();

  void preProcess() const {}
  void priorFilter(const GeoVaLs &) const {}
  void postFilter(const ioda::ObsVector &, const ObsDiagnostics &);

  const oops::Variables & requiredGeoVaLs() const {return nogeovals_;}
  const oops::Variables & requiredHdiagnostics() const {return extradiagvars_;}

 private:
  void print(std::ostream &) const;
  const eckit::LocalConfiguration config_;
  const oops::Variables nogeovals_;
  oops::Variables extradiagvars_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_OBSDIAGNOSTICSWRITER_H_
