/*
 * (C) Copyright 2026 NOAA/NWS/OMD
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_MARINE_SEAICEFREEBOARD_OBSSEAICEFREEBOARDPARAMETERS_H_
#define UFO_OPERATORS_MARINE_SEAICEFREEBOARD_OBSSEAICEFREEBOARDPARAMETERS_H_

#include "oops/util/parameters/NumericConstraints.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/ObsOperatorParametersBase.h"

namespace ufo {

/// The altimeter type observing sea ice freeboard.
///   - Lidar altimeters (e.g. ICESat-2/ATLAS) see the air/snow interface and so measure
///     total freeboard (ice + snow above sea level).
///   - Ka-band radar altimeters (e.g. SARAL/AltiKa) penetrate little into the snow and are
///     also assumed to see the air/snow interface (total freeboard).
///   - Ku-band radar altimeters (e.g. CryoSat-2/SIRAL) are assumed to penetrate to the
///     snow/ice interface and so measure ice-only freeboard.
/// The instruments differ only in the snow-thickness coefficient of the freeboard
/// equation -- see ObsSeaIceFreeboard.cc.
enum class SeaIceFreeboardInstrument {
  Lidar,
  RadarKa,
  RadarKu
};
struct SeaIceFreeboardInstrumentParameterTraitsHelper {
  typedef SeaIceFreeboardInstrument EnumType;
  static constexpr char enumTypeName[] = "SeaIceFreeboardInstrument";
  static constexpr util::NamedEnumerator<SeaIceFreeboardInstrument> namedValues[] = {
    { SeaIceFreeboardInstrument::Lidar, "lidar" },
    { SeaIceFreeboardInstrument::RadarKa, "radar Ka" },
    { SeaIceFreeboardInstrument::RadarKu, "radar Ku" }
  };
};

}  // namespace ufo

namespace oops {

template <>
struct ParameterTraits<ufo::SeaIceFreeboardInstrument> :
    public EnumParameterTraits<ufo::SeaIceFreeboardInstrumentParameterTraitsHelper>
{};

}  // namespace oops

namespace ufo {

// -----------------------------------------------------------------------------
/// SeaIceFreeboard observation operator Parameters class.
///
/// All densities are required (no hardcoded defaults) so every configuration is
/// explicit about the values used to convert model ice/snow thickness into freeboard.
/// The bounds below are realistic physical ranges; they also guarantee that sea water
/// is denser than both sea ice and snow, so the freeboard coefficients are well defined.
class ObsSeaIceFreeboardParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsSeaIceFreeboardParameters, ObsOperatorParametersBase)

 public:
  /// Altimeter type: "lidar" or "radar Ka" (total freeboard), or "radar Ku" (ice-only
  /// freeboard).
  oops::RequiredParameter<SeaIceFreeboardInstrument> instrument{"instrument", this};

  oops::RequiredParameter<double> seaIceDensity{
      "sea ice density", this, {oops::minConstraint(850.0), oops::maxConstraint(950.0)}};
  oops::RequiredParameter<double> surfaceSnowDensity{
      "surface snow density", this, {oops::minConstraint(50.0), oops::maxConstraint(450.0)}};
  oops::RequiredParameter<double> seaWaterDensity{
      "sea water density", this, {oops::minConstraint(990.0), oops::maxConstraint(1035.0)}};
};

}  // namespace ufo

#endif  // UFO_OPERATORS_MARINE_SEAICEFREEBOARD_OBSSEAICEFREEBOARDPARAMETERS_H_
