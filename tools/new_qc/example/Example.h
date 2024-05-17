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
#include "oops/base/ObsVariables.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
// TODO: modify the list of Parameter classes to include
// depending on what is used below.
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/interface/ObsFilterBase.h"
#include "tools/new_qc/example/Example.interface.h"
#include "ufo/filters/FilterParametersBase.h"
#include "ufo/ObsTraits.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

/// \brief Options controlling the operation of the Example filter.

class ExampleParameters : public oops::ObsFilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(ExampleParameters, oops::ObsFilterParametersBase)
  // TODO: list all parameters here.
  oops::Parameter<int> myParameter
    {"my Parameter",
     "This is a Parameter",
     4,
     this};

  oops::OptionalParameter<float> myOptionalParameter
    {"my OptionalParameter",
     "This is an OptionalParameter",
     this};

  oops::RequiredParameter<std::string> myRequiredParameter
    {"my RequiredParameter",
     "This is a RequiredParameter",
     this};

  oops::OptionalParameter<Variable> myVariableParameter
    {"my VariableParameter",
     "This is an OptionalParameter holding a Variable",
     this};
};

/// Example filter

class Example : public oops::interface::ObsFilterBase<ObsTraits>,
                private util::ObjectCounter<Example> {
 public:
  typedef ExampleParameters Parameters_;

  static const std::string classname() {return "ufo::Example";}

  Example(ioda::ObsSpace &, const Parameters_ &,
          std::shared_ptr<ioda::ObsDataVector<int> >,
          std::shared_ptr<ioda::ObsDataVector<float> >);
  ~Example();

  void preProcess() override {}
  void priorFilter(const GeoVaLs &) override;
  void postFilter(const GeoVaLs &,
                  const ioda::ObsVector &,
                  const ioda::ObsVector &,
                  const ObsDiagnostics &) override;
  void checkFilterData(const oops::FilterStage filterStage) override {};

  oops::Variables requiredVars() const override {return geovars_;}
  oops::ObsVariables requiredHdiagnostics() const override {return diagnostics_;}

 private:
  void print(std::ostream &) const override;
  F90check key_;

  ioda::ObsSpace & obsdb_;
  oops::Variables geovars_;
  oops::ObsVariables diagnostics_;
  ioda::ObsDataVector<int> & flags_;
};

}  // namespace ufo

#endif  // TOOLS_NEW_QC_EXAMPLE_EXAMPLE_H_
