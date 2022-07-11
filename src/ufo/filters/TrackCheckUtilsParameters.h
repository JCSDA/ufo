/*
 * (C) 2021 Crown Copyright Met Office. All rights reserved.
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_TRACKCHECKUTILSPARAMETERS_H_
#define UFO_FILTERS_TRACKCHECKUTILSPARAMETERS_H_

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "ufo/filters/FilterParametersBase.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

  enum class SurfaceObservationSubtype {
    LNDSYN, SHPSYN, BUOY, MOBSYN, OPENROAD, TEMP, BATHY, TESAC, BUOYPROF, LNDSYB, SHPSYB
  };
  struct SurfaceObservationSubtypeParameterTraitsHelper {
    typedef SurfaceObservationSubtype EnumType;
    static constexpr char enumTypeName[] = "SurfaceObservationSubtype";
    static constexpr util::NamedEnumerator<SurfaceObservationSubtype> namedValues[] = {
      { SurfaceObservationSubtype::LNDSYN, "LNDSYN" },
      { SurfaceObservationSubtype::SHPSYN, "SHPSYN" },
      { SurfaceObservationSubtype::BUOY, "BUOY" },
      { SurfaceObservationSubtype::MOBSYN, "MOBSYN" },
      { SurfaceObservationSubtype::OPENROAD, "OPENROAD" },
      { SurfaceObservationSubtype::TEMP, "TEMP" },
      { SurfaceObservationSubtype::BATHY, "BATHY" },
      { SurfaceObservationSubtype::TESAC, "TESAC" },
      { SurfaceObservationSubtype::BUOYPROF, "BUOYPROF" },
      { SurfaceObservationSubtype::LNDSYB, "LNDSYB" },
      { SurfaceObservationSubtype::SHPSYB, "SHPSYB" }
    };
  };
}  // namespace ufo

namespace oops {

template<>
struct ParameterTraits<ufo::SurfaceObservationSubtype> :
    public EnumParameterTraits<ufo::SurfaceObservationSubtypeParameterTraitsHelper>
{};

}  // namespace oops

namespace ufo {

/// \brief Options controlling the operation of the track check filter.
class TrackCheckUtilsParameters : public FilterParametersBase {
  OOPS_ABSTRACT_PARAMETERS(TrackCheckUtilsParameters, FilterParametersBase)

 public:
  /// Variable storing integer-valued or string-valued station IDs.
  /// Observations taken by each station are checked separately.
  ///
  /// If not set and observations were grouped into records when the observation space was
  /// constructed, each record is assumed to consist of observations taken by a separate
  /// station. If not set and observations were not grouped into records, all observations are
  /// assumed to have been taken by a single station.
  ///
  /// Note: the variable used to group observations into records can be set with the
  /// \c obs space.obsdatain.obsgrouping.groupvariable YAML option.
  oops::OptionalParameter<Variable> stationIdVariable{
    "station_id_variable", this};
};

}  // namespace ufo

#endif  // UFO_FILTERS_TRACKCHECKUTILSPARAMETERS_H_
