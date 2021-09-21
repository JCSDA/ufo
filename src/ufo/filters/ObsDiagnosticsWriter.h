/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_OBSDIAGNOSTICSWRITER_H_
#define UFO_FILTERS_OBSDIAGNOSTICSWRITER_H_

#include <memory>
#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/interface/ObsFilterBase.h"
#include "oops/util/Printable.h"
#include "ufo/filters/Variables.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/ObsTraits.h"

namespace eckit {
  class Configuration;
  class LocalConfiguration;
}

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
class GeoVaLs;

class ObsDiagnosticsWriter : public oops::interface::ObsFilterBase<ObsTraits> {
 public:
  ObsDiagnosticsWriter(ioda::ObsSpace &, const eckit::Configuration &,
                       std::shared_ptr<ioda::ObsDataVector<int> >,
                       std::shared_ptr<ioda::ObsDataVector<float> >);
  ~ObsDiagnosticsWriter() {}

  void preProcess() override {}
  void priorFilter(const GeoVaLs &) override {}
  void postFilter(const ioda::ObsVector &, const ioda::ObsVector &,
                  const ObsDiagnostics & diags) override {
    diags.write(config_);
  }

  oops::Variables requiredVars() const override {return nogeovals_;}
  oops::Variables requiredHdiagnostics() const override {return extradiagvars_;}

 private:
  void print(std::ostream &) const override;
  const eckit::LocalConfiguration config_;
  const oops::Variables nogeovals_;
  oops::Variables extradiagvars_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_OBSDIAGNOSTICSWRITER_H_
