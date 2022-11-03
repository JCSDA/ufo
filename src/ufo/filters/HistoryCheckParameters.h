/*
 * (C) Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#ifndef UFO_FILTERS_HISTORYCHECKPARAMETERS_H_
#define UFO_FILTERS_HISTORYCHECKPARAMETERS_H_

#include <utility>

#include "ioda/ObsSpaceParameters.h"

#include "oops/util/Duration.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/ParameterTraits.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/StuckCheckParameters.h"
#include "ufo/filters/TrackCheckShipParameters.h"
#include "ufo/filters/TrackCheckUtilsParameters.h"

namespace ufo {

/// \brief Options controlling the operation of history check filter.
class HistoryCheckParameters : public TrackCheckUtilsParameters {
  OOPS_CONCRETE_PARAMETERS(HistoryCheckParameters, TrackCheckUtilsParameters)

 public:
    /// Surface observation subtype determining if the track and/or stuck check could be run
    /// \todo It is possible to use this to decide on the window length. Should the manual entry
    /// method be kept?
    oops::RequiredParameter<SurfaceObservationSubtype> surfaceObservationSubtype{
      "input category", this};

    /// Amount of time before start of assimilation window to collect for the history check
    oops::RequiredParameter<util::Duration> timeBeforeStartOfWindow {
      "time before start of window", this
    };

    /// Amount of time (default zero) after end of assimilation window to collect
    /// for the history check
    oops::Parameter<util::Duration> timeAfterEndOfWindow {
      "time after end of window", util::Duration("PT0S"), this
    };

    /// The options for running the ship track check filter which can be optionally run, should the
    /// subtype not be LNDSYN/LNDSYB. These must be filled in in order for the track filter to run.
    oops::OptionalParameter<TrackCheckShipCoreParameters> trackCheckShipParameters {
      "ship track check parameters", this
    };

    /// The options for running the stuck check filter which can be optionally run, should the
    /// subtype not be TEMP/BATHY/TESAC/BUOYPROF. These must be filled in in order for the stuck
    /// filter to run.
    oops::OptionalParameter<StuckCheckCoreParameters> stuckCheckParameters {
      "stuck check parameters", this
    };

    /// Creates a new obs space with the wider window that is determined by the observation subtype.
    /// Needs: name (can be set with setValue), simulated variables, obsdatain.obsfile.
    oops::RequiredParameter<ioda::ObsTopLevelParameters> largerObsSpace {
      "obs space", this
    };

    /// Controls whether all of the larger obs space's variables are reset to match the primary
    /// obs space's when filter is run. Used for unit testing (esp. stuck check portions).
    oops::Parameter<bool> resetLargerObsSpaceVariables {
      "reset larger obs space variables", false, this
    };

    /// Maximum number of characters for a string-labelled station id.
    /// This is used to ensure unique integer hashes if station ids are string labels.
    oops::Parameter<int> stationIdMaxStringLength {
      "station id max string length", 24, this
    };
};

}  // namespace ufo

#endif  // UFO_FILTERS_HISTORYCHECKPARAMETERS_H_
