/*
 * (C) Copyright 2025 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/PotentialTemperatureFromTemperature.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/Constants.h"

namespace ufo {

static ObsFunctionMaker<PotentialTemperatureFromTemperature>
  makerPotentialTemperatureFromTemperature_("PotentialTemperatureFromTemperature");

PotentialTemperatureFromTemperature::PotentialTemperatureFromTemperature(
  const eckit::LocalConfiguration & conf
): invars_() {
  oops::Log::trace() << "PotentialTemperatureFromTemperature constructor" << std::endl;
  // Initialize options
  options_.deserialize(conf);
  std::ostringstream errString;
  if (options_.useSurfacePressure.value() == false) {
    ASSERT(options_.geovalsPressure.value() != boost::none);
  }

  if (options_.useSurfacePressure.value()) {
     invars_ += Variable("GeoVaLs/air_pressure_at_surface");
     invars_ += Variable("GeoVaLs/average_surface_temperature_within_field_of_view");
  } else {
     invars_ += Variable("GeoVaLs/air_pressure");
     invars_ += Variable("GeoVaLs/air_potential_temperature");
  }
}

// -----------------------------------------------------------------------------

PotentialTemperatureFromTemperature::~PotentialTemperatureFromTemperature() {
  oops::Log::trace() << "PotentialTemperatureFromTemperature destructor" << std::endl;
}

// -----------------------------------------------------------------------------

void PotentialTemperatureFromTemperature::compute(const ObsFilterData & in,
                                    ioda::ObsDataVector<float> & out) const {
  oops::Log::trace() << "PotentialTemperatureFromTemperature compute start" << std::endl;
  // Get dimension
  const size_t nlocs = in.nlocs();
  std::vector<float> potentialTemperature(nlocs, 0.0);
  std::vector<float> pressure(nlocs, 1.0E+5);
  std::vector<float> temperature(nlocs, 0.0);

  if (options_.useSurfacePressure.value()) {
     in.get(Variable("GeoVaLs/air_pressure_at_surface"), pressure);
     in.get(Variable("GeoVaLs/average_surface_temperature_within_field_of_view"), temperature);
     // Calculate potential temperature
     for (size_t iloc = 0; iloc < nlocs; ++iloc) {
        potentialTemperature[iloc] =  temperature[iloc] * pow(Constants::pref / pressure[iloc],
                     Constants::rd_over_cp);
     }
  } else {
     float pressure_value = options_.geovalsPressure.value().get();
     const size_t nlevs = in.nlevs(Variable("GeoVaLs/air_pressure"));
     std::vector<float> dprs_min(nlocs, 9999.0);
     std::vector<float> pre_levl(nlocs);
     std::vector<float> air_potential_temperature(nlocs);
     for (size_t ilev = 0; ilev < nlevs; ++ilev) {
        in.get(Variable("GeoVaLs/air_pressure"), ilev, pre_levl);
        in.get(Variable("GeoVaLs/air_potential_temperature"), ilev, air_potential_temperature);
        for (size_t iloc = 0; iloc < nlocs; ++iloc) {
           float pressure_diff = std::fabs(pre_levl[iloc] - pressure_value);
           if (pressure_diff < dprs_min[iloc]) {
              dprs_min[iloc] = pressure_diff;
              potentialTemperature[iloc] = air_potential_temperature[iloc];
           }
        }
     }
  }
  out[0] = potentialTemperature;
  oops::Log::trace() << "PotentialTemperatureFromTemperature compute complete" << std::endl;
}

// -----------------------------------------------------------------------------

const ufo::Variables & PotentialTemperatureFromTemperature::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
