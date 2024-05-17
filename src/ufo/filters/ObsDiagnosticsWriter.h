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
#include <vector>

#include "oops/base/ObsVariables.h"
#include "oops/base/Variables.h"
#include "oops/interface/ObsFilterBase.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/ObsTraits.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
  class ObsVector;
}

namespace ufo {

/// \brief Parameters controlling ObsDiagnosticsWriter
class ObsDiagnosticsWriterParameters : public oops::ObsFilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsDiagnosticsWriterParameters, oops::ObsFilterParametersBase)
 public:
  oops::OptionalParameter<std::vector<Variable>> filterVariables{"filter variables", this};

  GeoVaLsParameters diags{this};
};


class ObsDiagnosticsWriter : public oops::interface::ObsFilterBase<ObsTraits> {
 public:
  typedef ObsDiagnosticsWriterParameters Parameters_;
  ObsDiagnosticsWriter(ioda::ObsSpace &, const Parameters_ &,
                       std::shared_ptr<ioda::ObsDataVector<int> >,
                       std::shared_ptr<ioda::ObsDataVector<float> >);
  ~ObsDiagnosticsWriter() {}

  void preProcess() override {}
  void priorFilter(const GeoVaLs &) override {}
  void postFilter(const GeoVaLs &,
                  const ioda::ObsVector &,
                  const ioda::ObsVector &,
                  const ObsDiagnostics & diags) override {
    diags.write(params_.diags.toConfiguration());
  }
  void checkFilterData(const oops::FilterStage filterStage) override {}

  oops::Variables requiredVars() const override {return nogeovals_;}
  oops::ObsVariables requiredHdiagnostics() const override {return extradiagvars_;}

 private:
  void print(std::ostream &) const override;
  Parameters_ params_;
  const oops::Variables nogeovals_;
  oops::ObsVariables extradiagvars_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_OBSDIAGNOSTICSWRITER_H_
