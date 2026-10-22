/*
 * (C) Copyright 2026 NOAA/NWS/OMD
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_MARINE_SEAICEFREEBOARD_OBSSEAICEFREEBOARD_H_
#define UFO_OPERATORS_MARINE_SEAICEFREEBOARD_OBSSEAICEFREEBOARD_H_

#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/ObsOperatorBase.h"
#include "ufo/operators/marine/seaicefreeboard/ObsSeaIceFreeboardParameters.h"

/// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// Setup and checks shared by the SeaIceFreeboard nonlinear and linear operators.
namespace seaicefreeboard {

/// Coefficients of the freeboard equation F = alpha * h_ice + beta * h_snow.
struct Coefficients {
  double alpha;
  double beta;
};

/// Computes the freeboard coefficients from the densities and the instrument.
Coefficients coefficients(const ObsSeaIceFreeboardParameters &);

/// The geovals required by the operator: sea ice and snow thickness.
oops::Variables requiredVars();

/// Throws unless "seaIceFreeboard" is the only simulated variable.
void checkSimulatedVariable(const ioda::ObsSpace &);

/// Throws unless the ice and snow thickness geovals are single-category.
void checkSingleCategory(const GeoVaLs &);

}  // namespace seaicefreeboard

// -----------------------------------------------------------------------------
/// Sea ice freeboard observation operator, covering both lidar and radar altimetry.
///
/// Computes freeboard F = alpha * h_ice + beta * h_snow, where
/// alpha = (rho_w - rho_i) / rho_w. The `instrument` Parameter selects beta:
///   - lidar, radar Ka (total freeboard, air/snow interface):
///     beta = (rho_w - rho_s) / rho_w
///   - radar Ku (ice-only freeboard, snow/ice interface): beta = -rho_s / rho_w
///     The slower wave speed in snow is ignored: the observation must be the ice
///     freeboard already corrected for it, not the raw radar freeboard.
///
/// Freeboard is an intensive property of the ice surface: a height measured over
/// ice, not a grid-cell mean. The geovals `sea_ice_thickness` and
/// `sea_ice_snow_thickness` must be single-category thicknesses per unit *ice*
/// area, so no sea ice concentration weighting is applied and
/// `sea_ice_area_fraction` is not required. Screen out low-concentration points
/// with a Domain Check filter on GeoVaLs/sea_ice_area_fraction, where freeboard
/// is undefined.
class ObsSeaIceFreeboard : public ObsOperatorBase,
                          private util::ObjectCounter<ObsSeaIceFreeboard> {
 public:
  typedef ObsSeaIceFreeboardParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;

  static const std::string classname() {return "ufo::ObsSeaIceFreeboard";}

  ObsSeaIceFreeboard(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsSeaIceFreeboard();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t &) const override;

// Other
  const oops::Variables & requiredVars() const override {return varin_;}

 private:
  void print(std::ostream &) const override;
  const oops::Variables varin_;
  const seaicefreeboard::Coefficients coefs_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_MARINE_SEAICEFREEBOARD_OBSSEAICEFREEBOARD_H_
