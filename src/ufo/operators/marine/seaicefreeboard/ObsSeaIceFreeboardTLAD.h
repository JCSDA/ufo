/*
 * (C) Copyright 2026 NOAA/NWS/OMD
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_MARINE_SEAICEFREEBOARD_OBSSEAICEFREEBOARDTLAD_H_
#define UFO_OPERATORS_MARINE_SEAICEFREEBOARD_OBSSEAICEFREEBOARDTLAD_H_

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/LinearObsOperatorBase.h"
#include "ufo/operators/marine/seaicefreeboard/ObsSeaIceFreeboard.h"

// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// Sea ice freeboard TL/AD observation operator (lidar and radar altimetry).
///
/// F = alpha * h_ice + beta * h_snow is exactly linear in the thicknesses, so the
/// tangent linear is the same expression applied to the increments and no
/// trajectory values are needed - only which locations the trajectory marked
/// missing. Freeboard carries no sensitivity to sea ice concentration.
class ObsSeaIceFreeboardTLAD : public LinearObsOperatorBase,
                               private util::ObjectCounter<ObsSeaIceFreeboardTLAD> {
 public:
  typedef ObsSeaIceFreeboardParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;

  static const std::string classname() {return "ufo::ObsSeaIceFreeboardTLAD";}

  ObsSeaIceFreeboardTLAD(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsSeaIceFreeboardTLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, ObsDiagnostics &, const QCFlags_t &) override;
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &) const override;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &) const override;

  // Other
  const oops::Variables & requiredVars() const override {return varin_;}

 private:
  void print(std::ostream &) const override;
  const oops::Variables varin_;
  const seaicefreeboard::Coefficients coefs_;

  // Locations where the trajectory had a missing thickness, saved by setTrajectory.
  std::vector<bool> missingTraj_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_MARINE_SEAICEFREEBOARD_OBSSEAICEFREEBOARDTLAD_H_
