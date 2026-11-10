/*
 * (C) Copyright 2026 NOAA/NWS/OMD
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/marine/seaicefreeboard/ObsSeaIceFreeboard.h"

#include <ostream>
#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/GeoVaLs.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsSeaIceFreeboard> makerSeaIceFreeboard_("SeaIceFreeboard");
// -----------------------------------------------------------------------------

namespace seaicefreeboard {

Coefficients coefficients(const ObsSeaIceFreeboardParameters & params) {
  const double rhoIce = params.seaIceDensity.value();
  const double rhoSnow = params.surfaceSnowDensity.value();
  const double rhoWater = params.seaWaterDensity.value();

  Coefficients coefs;
  coefs.alpha = (rhoWater - rhoIce) / rhoWater;
  switch (params.instrument.value()) {
    case SeaIceFreeboardInstrument::Lidar:
    case SeaIceFreeboardInstrument::RadarKa:
      // Total freeboard: air/snow interface above sea level.
      coefs.beta = (rhoWater - rhoSnow) / rhoWater;
      break;
    case SeaIceFreeboardInstrument::RadarKu:
      // Ice-only freeboard: snow/ice interface above sea level.
      // The slower Ku-band wave speed in snow is ignored, so the observation must be
      // the ice freeboard already corrected for it, not the raw radar freeboard. The
      // raw radar freeboard would use
      //   beta = -((c/c_s - 1) * rho_w + rho_s) / rho_w,  c/c_s = (1 + 0.00051 rho_s)^1.5
      coefs.beta = -rhoSnow / rhoWater;
      break;
    default:
      throw eckit::BadValue("SeaIceFreeboard: unsupported instrument", Here());
  }
  return coefs;
}

// -----------------------------------------------------------------------------

oops::Variables requiredVars() {
  return oops::Variables(std::vector<std::string>{"sea_ice_thickness",
                                                  "sea_ice_snow_thickness"});
}

// -----------------------------------------------------------------------------

void checkSimulatedVariable(const ioda::ObsSpace & odb) {
  if (!(odb.assimvariables().size() == 1 &&
        odb.assimvariables().variables()[0] == "seaIceFreeboard")) {
    throw eckit::BadValue(
      "SeaIceFreeboard can only simulate the variable \"seaIceFreeboard\"", Here());
  }
}

// -----------------------------------------------------------------------------

void checkSingleCategory(const GeoVaLs & gv) {
  for (const char * var : {"sea_ice_thickness", "sea_ice_snow_thickness"}) {
    const size_t ncat = gv.nlevs(oops::Variable{var});
    if (ncat != 1) {
      throw eckit::BadValue(
        std::string("SeaIceFreeboard expects single-category \"") + var + "\" (per unit ice "
        "area), but the geovals have " + std::to_string(ncat) + " categories", Here());
    }
  }
}

}  // namespace seaicefreeboard

// -----------------------------------------------------------------------------

ObsSeaIceFreeboard::ObsSeaIceFreeboard(const ioda::ObsSpace & odb, const Parameters_ & params)
  : ObsOperatorBase(odb), varin_(seaicefreeboard::requiredVars()),
    coefs_(seaicefreeboard::coefficients(params))
{
  seaicefreeboard::checkSimulatedVariable(odb);
  oops::Log::trace() << "ObsSeaIceFreeboard constructor done." << std::endl;
}

// -----------------------------------------------------------------------------

ObsSeaIceFreeboard::~ObsSeaIceFreeboard() {
  oops::Log::trace() << "ObsSeaIceFreeboard destructor done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaIceFreeboard::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                                     ObsDiagnostics &, const QCFlags_t &) const {
  const double missing = util::missingValue<double>();
  const size_t nlocs = ovec.size();

  seaicefreeboard::checkSingleCategory(gv);

  std::vector<double> iceThick(nlocs);
  std::vector<double> snowThick(nlocs);
  gv.getAtLevel(iceThick, oops::Variable{"sea_ice_thickness"}, 0);
  gv.getAtLevel(snowThick, oops::Variable{"sea_ice_snow_thickness"}, 0);

  for (size_t i = 0; i < nlocs; ++i) {
    if (iceThick[i] == missing || snowThick[i] == missing) {
      ovec[i] = missing;
    } else {
      ovec[i] = coefs_.alpha * iceThick[i] + coefs_.beta * snowThick[i];
    }
  }

  oops::Log::trace() << "ObsSeaIceFreeboard::simulateObs done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSeaIceFreeboard::print(std::ostream & os) const {
  os << "ObsSeaIceFreeboard operator, alpha=" << coefs_.alpha << ", beta=" << coefs_.beta;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
