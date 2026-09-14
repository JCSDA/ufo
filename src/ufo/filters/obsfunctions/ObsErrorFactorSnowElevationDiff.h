/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSERRORFACTORSNOWELEVATIONDIFF_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSERRORFACTORSNOWELEVATIONDIFF_H_

#include <memory>
#include <string>
#include <vector>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

/// \brief Options controlling ObsErrorFactorSnowElevationDiff ObsFunction
class ObsErrorFactorSnowElevationDiffParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsErrorFactorSnowElevationDiffParameters, Parameters)

 public:
  oops::RequiredParameter<std::string> obs_elevation_var{"obs_elevation_var", this};
  oops::RequiredParameter<std::string> model_elevation_var{"model_elevation_var", this};
  oops::RequiredParameter<float> elevation_scale_h{"elevation_scale_h", this};
};

// -----------------------------------------------------------------------------

/// \brief Inflate the observation error based on elevation difference between model and observation.
///
/// This routine computes an observation error inflation factor based on the elevation difference
/// between the model surface elevation and the observed station elevation.
/// The inflation factor is computed as: 1 / exp(-1 * dz^2 / (h^2))
/// where dz = |model_elevation - obs_elevation| and h is the elevation_scale_h parameter.
///
/// ~~~~
///
/// ### example configurations for a FilterBase derived class: ###
///
///     - filter: BlackList
///       filter variables:
///       - name: snowDepth
///       action:
///         name: inflate error
///         inflation variable:
///           name: ObsFunction/ObsErrorFactorSnowElevationDiff
///           options:
///             obs_elevation_var: MetaData/stationElevation
///             model_elevation_var: GeoVaLs/filtered_orography
///             elevation_scale_h: 100.0
///
class ObsErrorFactorSnowElevationDiff : public ObsFunctionBase<float> {
 public:
  static const std::string classname() {return "ObsErrorFactorSnowElevationDiff";}

  explicit ObsErrorFactorSnowElevationDiff(const eckit::Configuration &config);
  ~ObsErrorFactorSnowElevationDiff();

  void compute(const ObsFilterData &, ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  std::unique_ptr<ObsErrorFactorSnowElevationDiffParameters> options_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSERRORFACTORSNOWELEVATIONDIFF_H_
