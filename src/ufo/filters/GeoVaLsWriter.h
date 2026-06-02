/*
 * (C) Copyright 2017-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/ObsVariables.h"
#include "oops/base/Variables.h"
#include "ufo/filters/GenericFilterParameters.h"
#include "ufo/ObsFilterBase.h"

namespace ioda {
class ObsSpace;
class ObsVector;
template <typename DATA> class ObsDataVector;
}

namespace ufo {
class GeoVaLs;
class ObsDiagnostics;

// -----------------------------------------------------------------------------

class GeoVaLsWriter : public ObsFilterBase {
  template <typename DATA> using ObsDataPtr_ = std::shared_ptr<ioda::ObsDataVector<DATA> >;

 public:
  typedef GenericFilterParameters       Parameters_;

  GeoVaLsWriter(const ioda::ObsSpace &, const Parameters_ &,
                ObsDataPtr_<int>, ObsDataPtr_<float>);
  ~GeoVaLsWriter() = default;

  void preProcess() override {}

  void priorFilter(const GeoVaLs & gv) override;

  void postFilter(const GeoVaLs &,
                  const ioda::ObsVector &,
                  const ioda::ObsVector &,
                  const ObsDiagnostics &) override;

  void checkFilterData(const FilterStage filterStage) override;

  oops::Variables requiredVars() const override {return novars_;};
  oops::ObsVariables requiredHdiagnostics() const override {return noobsvars_;};

 private:
  const eckit::LocalConfiguration config_;
  const oops::Variables novars_;  // could be used to determine what needs saving
  const oops::ObsVariables noobsvars_;
  void print(std::ostream &) const override;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
