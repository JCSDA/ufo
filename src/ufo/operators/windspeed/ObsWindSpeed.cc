/*
 * (C) Crown Copyright 2024 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/operators/windspeed/ObsWindSpeed.h"

#include <cmath>
#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/ObsOperatorBase.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsWindSpeed> makerWindSpeed_("WindSpeed");
// -----------------------------------------------------------------------------

ObsWindSpeed::ObsWindSpeed(const ioda::ObsSpace & odb,
                 const Parameters_ &params)
  : ObsOperatorBase(odb), model_surface_eastward_wind_{params.model_surface_eastward_wind.value()},
    model_surface_northward_wind_{params.model_surface_northward_wind.value()},
    varin_{{model_surface_eastward_wind_, model_surface_northward_wind_}}
{
  oops::Log::trace() << "ObsWindSpeed created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsWindSpeed::~ObsWindSpeed() {
  oops::Log::trace() << "ObsWindSpeed destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsWindSpeed::simulateObs(const GeoVaLs & geovals, ioda::ObsVector & hofx,
                             ObsDiagnostics &, const QCFlags_t&) const {
// Observation operator

  // Check hofx size and initialise hofx to zero:
  ASSERT(geovals.nlocs() == hofx.nlocs());
  hofx.zero();
  // Get number of locations
  std::size_t nwinds = geovals.nlocs();

  // Get u and v winds

  std::vector<double> u_wind(nwinds);  // u wind component
  geovals.get(u_wind, model_surface_eastward_wind_);
  std::vector<double> v_wind(nwinds);  // v wind component
  geovals.get(v_wind, model_surface_northward_wind_);

  // calculate wind speed
  for (size_t w = 0; w < nwinds; w++) {
    hofx[w] += sqrt(pow(u_wind[w], 2.0) + pow(v_wind[w], 2.0));
  }

  oops::Log::trace() << "ObsWindSpeed: observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsWindSpeed::print(std::ostream & os) const {
  os << "ObsWindSpeed::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
