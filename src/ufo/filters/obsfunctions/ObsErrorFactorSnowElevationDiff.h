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
  oops::Parameter<std::string> obs_elevation_var{"observation_elevation", "MetaData/stationElevation", this};
  oops::Parameter<std::string> model_elevation_var{"model_elevation", "GeoVaLs/filtered_orography", this};
  oops::Parameter<float> elevation_scale_m{"elevation_scale_m", 800., this};
};

// -----------------------------------------------------------------------------

/// \brief Inflate observation error based on elevation difference between model and observation.
///
/// This routine computes an observation error inflation factor based on the elevation difference
/// between the model surface elevation and the observed station elevation.
/// The inflation factor is computed as: 1 / exp(-1 * dz^2 / (h^2))
/// where dz = |model_elevation - obs_elevation| 
/// and h is the elevation_scale parameter in m (Elevation difference at which the obs-error will be inflated by a factor of exp(1)).
///
/// Authors
/// First draft: Github Copilot. 
/// Revisiewd and tested, in accordance with NOAA's use of AI tools, by Tseganeh Z. Gichamo
/// 
/// ### example configurations for application of this filter: ###
///
///     - filter: Perform Action
///       filter variables:
///       - name: totalSnowDepth
///       action:
///         name: inflate error
///         inflation variable:
///           name: ObsFunction/ObsErrorFactorSnowElevationDiff
///           options:
///             observation_elevation: MetaData/stationElevation
///             model_elevation: GeoVaLs/filtered_orography
///             elevation_scale_m: 800.0   // (Unit m) 
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
