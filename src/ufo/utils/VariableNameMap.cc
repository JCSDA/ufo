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
  Aliases_.emplace("airTemperature", oops::Variable{"air_temperature"});
  Aliases_.emplace("windEastward", oops::Variable{"eastward_wind"});
  Aliases_.emplace("windNorthward", oops::Variable{"northward_wind"});
  Aliases_.emplace("specificHumidity", oops::Variable{"specific_humidity"});
  Aliases_.emplace("relativeHumidity", oops::Variable{"relative_humidity"});
  Aliases_.emplace("pressure", oops::Variable{"air_pressure"});
  Aliases_.emplace("virtualTemperature", oops::Variable{"virtual_temperature"});
  Aliases_.emplace("potentialTemperature", oops::Variable{"potential_temperature"});
  Aliases_.emplace("stationPressure", oops::Variable{"surface_pressure"});
  Aliases_.emplace("depthBelowWaterSurface", oops::Variable{"ocean_depth"});
  Aliases_.emplace("waterTemperature", oops::Variable{"ocean_temperature"});
  Aliases_.emplace("surfacePressure", oops::Variable{"surface_pressure"});
  Aliases_.emplace("seaSurfaceTemperature", oops::Variable{"sea_surface_temperature"});
  Aliases_.emplace("equivalentReflectivityFactor", oops::Variable
                                                                {"equivalent_reflectivity_factor"});
  Aliases_.emplace("salinity", oops::Variable{"sea_water_salinity"});
  Aliases_.emplace("seaSurfaceSalinity", oops::Variable{"sea_surface_salinity"});
  Aliases_.emplace("ozoneProfile", oops::Variable{"mole_fraction_of_ozone_in_air"});
  Aliases_.emplace("seaIceFraction", oops::Variable{"sea_ice_area_fraction"});

// The lines below (and yaml files) should be removed as all name are hard-wired because of
// the obs and model naming conventions
  if (aliasFile.is_initialized()) {
    eckit::YAMLConfiguration conf(eckit::PathName(aliasFile.value()));
    VariableMapParameters mappingParams;
    mappingParams.validateAndDeserialize(conf);

    for (VariableNameParameters const& variableMap : *mappingParams.variableMap.value()) {
      Aliases_.insert_or_assign(variableMap.name, oops::Variable{variableMap.alias});
    }
  }
}

VariableNameMap::~VariableNameMap() {}

oops::Variable VariableNameMap::convertName(const std::string & name) const {
  const auto & it = Aliases_.find(name);
  if (it == Aliases_.end()) {
    return oops::Variable{name};
  } else {
    return it->second;
  }
}

oops::Variables VariableNameMap::convertName(const oops::ObsVariables & vars) const {
  oops::Variables aliasvars;
  for (size_t jv = 0; jv < vars.size(); ++jv) {
    aliasvars.push_back(convertName(vars[jv]));
  }
  return aliasvars;
}
}  // namespace ufo
