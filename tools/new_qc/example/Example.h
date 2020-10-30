/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef TOOLS_NEW_QC_EXAMPLE_EXAMPLE_H_
#define TOOLS_NEW_QC_EXAMPLE_EXAMPLE_H_

#include <memory>
#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "tools/new_qc/example/Example.interface.h"

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

/// Example filter

class Example : public util::Printable,
                private util::ObjectCounter<Example> {
 public:
  static const std::string classname() {return "ufo::Example";}

  Example(ioda::ObsSpace &, const eckit::Configuration &,
          std::shared_ptr<ioda::ObsDataVector<int> >,
          std::shared_ptr<ioda::ObsDataVector<float> >);
  ~Example();

  void preProcess() const {}
  void priorFilter(const GeoVaLs &) const;
  void postFilter(const ioda::ObsVector &, const ObsDiagnostics &) const;

  const oops::Variables & requiredVars() const {return geovars_;}
  const oops::Variables & requiredHdiagnostics() const {return diagnostics_;}

 private:
  void print(std::ostream &) const;
  F90check key_;

  ioda::ObsSpace & obsdb_;
  oops::Variables geovars_;
  oops::Variables diagnostics_;
  ioda::ObsDataVector<int> & flags_;
};

}  // namespace ufo

#endif  // TOOLS_NEW_QC_EXAMPLE_EXAMPLE_H_
