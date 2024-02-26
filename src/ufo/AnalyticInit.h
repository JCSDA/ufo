/*
 * (C) Copyright 2020-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_ANALYTICINIT_H_
#define UFO_ANALYTICINIT_H_

#include <string>

#include "oops/interface/AnalyticInitBase.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/ObsTraits.h"

namespace ufo {
  class GeoVaLs;
  class SampledLocations;

/// Parameters for Analytic init (empty except for analytic init method defined
/// in the base class)
class AnalyticInitParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(AnalyticInitParameters, Parameters)
  oops::RequiredParameter<std::string> method{"method", this};
};

/// AnalyticInit: filling GeoVaLs with analytic formula
class AnalyticInit : public oops::interface::AnalyticInitBase<ObsTraits> {
 public:
  typedef AnalyticInitParameters Parameters_;
  explicit AnalyticInit(const eckit::Configuration &);
  void fillGeoVaLs(const SampledLocations &, GeoVaLs &) const override;

 private:
  Parameters_ options_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_ANALYTICINIT_H_
