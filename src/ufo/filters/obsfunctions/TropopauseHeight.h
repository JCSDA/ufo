/*
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_TROPOPAUSEHEIGHT_H_
#define UFO_FILTERS_OBSFUNCTIONS_TROPOPAUSEHEIGHT_H_

#include <string>
#include <vector>

#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"

namespace ufo {

/// \brief Options controlling the TropopauseHeight ObsFunction
class TropopauseHeightParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(TropopauseHeightParameters, Parameters)

 public:
  oops::Parameter<float> trop_max_pressure
    {"tropopause maximum pressure",
     100000.0f,
     this};

  oops::Parameter<float> trop_min_pressure
    {"tropopause minimum pressure",
     100.0f,
     this};

  oops::Parameter<float> trop_relative_humidity
    {"tropopause relative humidity",
     0.5f,
     this};

  oops::Parameter<bool> rhmethod
    {"rhmethod",
     "Use relative humidity method to determine tropopause height",
     false,
     this};
};

// -----------------------------------------------------------------------------

/// \brief Calculate the tropopause height in terms of model level.
///
/// The first method calculates the tropopause lapse rate in terms of log(pressure) coordinates
/// and compares that to the WMO definition of tropopause lapse rate. Once the calculated lapse
/// falls below the WMO definition, and provided the pressure at this level is between the
/// tropopause maximum and minimum pressures (these are user-defined), the tropopause level is set.
///
/// If the first method fails, the second method is called. This method first calculates the
/// relative humidity at each level. It then loops down through the atmosphere and at each level,
/// compares the relative humidity to the user-defined relative humidity at the tropopause. If
/// either relative humidity rises above the tropopause value, or if the pressure at that level
/// exceeds the maximum pressure of the tropopause, the tropopause level is set.
class TropopauseHeight : public ObsFunctionBase<int> {
 public:
  explicit TropopauseHeight(const eckit::LocalConfiguration &);
  ~TropopauseHeight();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<int> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  TropopauseHeightParameters options_;
  ufo::Variables invars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_TROPOPAUSEHEIGHT_H_
