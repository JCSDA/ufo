/*
 * (C) Copyright 2022 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/utils/VariableNameMap.h"
#include "eckit/config/YAMLConfiguration.h"
#include "eckit/filesystem/PathName.h"
#include "oops/util/Logger.h"

namespace ufo {
VariableNameMap::VariableNameMap(const boost::optional<std::string> & aliasFile) {
  Aliases_["airTemperature"] = "air_temperature";
  Aliases_["windEastward"] = "eastward_wind";
  Aliases_["windNorthward"] = "northward_wind";
  Aliases_["specificHumidity"] = "specific_humidity";
  Aliases_["relativeHumidity"] = "relative_humidity";
  Aliases_["pressure"] = "air_pressure";
  Aliases_["virtualTemperature"] = "virtual_temperature";
  Aliases_["potentialTemperature"] = "theta";
  Aliases_["stationPressure"] = "surface_pressure";
  Aliases_["depthBelowWaterSurface"] = "ocean_depth";
  Aliases_["waterTemperature"] = "ocean_temperature";
  Aliases_["surfacePressure"] = "surface_pressure";
  Aliases_["seaSurfaceTemperature"] = "sea_surface_temperature";
  Aliases_["equivalentReflectivityFactor"] = "equivalent_reflectivity_factor";
  Aliases_["salinity"] = "sea_water_salinity";
  Aliases_["seaSurfaceSalinity"] = "sea_surface_salinity";
  Aliases_["ozoneProfile"] = "mole_fraction_of_ozone_in_air";
  Aliases_["seaIceFraction"] = "sea_ice_area_fraction";

// The lines below (and yaml files) should be removed as all name are hard-wired because of
// the obs and model naming conventions
  if (aliasFile.is_initialized()) {
    eckit::YAMLConfiguration conf(eckit::PathName(aliasFile.value()));
    VariableMapParameters mappingParams;
    mappingParams.validateAndDeserialize(conf);

    for (VariableNameParameters const& variableMap : *mappingParams.variableMap.value()) {
      Aliases_[variableMap.name] = variableMap.alias;
    }
  }
}

VariableNameMap::~VariableNameMap() {}

std::string VariableNameMap::convertName(const std::string & name) const {
  std::string alias;
  const auto & it = Aliases_.find(name);
  if (it == Aliases_.end()) {
    return name;
  } else {
    alias = it->second;
  }
  return alias;
}

oops::Variables VariableNameMap::convertName(const oops::Variables & vars) const {
  std::vector<std::string> newvars;
  for (size_t jv = 0; jv < vars.size(); ++jv) {
    std::string alias = convertName(vars[jv]);
    newvars.push_back(alias);
  }
  oops::Variables aliasvars(newvars);
  return aliasvars;
}
}  // namespace ufo
