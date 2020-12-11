/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_METOFFICEBUDDYCHECKPARAMETERS_H_
#define UFO_FILTERS_METOFFICEBUDDYCHECKPARAMETERS_H_

#include <map>
#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "ufo/utils/Constants.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// \brief A box covering a specified (closed) interval of latitudes and longitudes.
class LatLonBoxParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(LatLonBoxParameters, Parameters)

 public:
  bool contains(float latitude, float longitude) const {
    return minLatitude <= latitude && latitude <= maxLatitude &&
           minLongitude <= longitude && longitude <= maxLongitude;
  }

  oops::Parameter<float> minLatitude{"min_latitude", -90, this};
  oops::Parameter<float> maxLatitude{"max_latitude", 90, this};
  oops::Parameter<float> minLongitude{"min_longitude", -180, this};
  oops::Parameter<float> maxLongitude{"max_longitude", 180, this};
};

/// \brief Options controlling the operation of the MetOfficeBuddyCheck filter.
class MetOfficeBuddyCheckParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(MetOfficeBuddyCheckParameters, Parameters)

 public:
  /// \name Parameters controlling buddy pair identification
  /// @{

  /// Maximum distance between two observations that may be classified as buddies, in km.
  oops::Parameter<float> searchRadius{"search_radius", 100, this};

  /// Variable storing string-valued or integer-valued station IDs.
  ///
  /// If not set and observations were grouped into records when the observation space was
  /// constructed, each record is assumed to consist of observations taken by a separate
  /// station. If not set and observations were not grouped into records, all observations are
  /// assumed to have been taken by a single station.
  ///
  /// Note: the variable used to group observations into records can be set with the
  /// \c obs space.obsdatain.obsgrouping.groupvariable YAML option.
  oops::OptionalParameter<Variable> stationIdVariable{"station_id_variable", this};

  /// Number of zonal bands to split the Earth's surface into when building a search data structure.
  ///
  /// Apart from the impact on the speed of buddy identification, this parameter also affects the
  /// order in which observations are processed and thus the final estimates of gross error
  /// probabilities, since the probability updates made when checking individual observation pairs
  /// are not commutative.
  oops::Parameter<int> numZonalBands{"num_zonal_bands", 24, this};

  /// Whether to include pressure in the sorting criteria used when building a search data
  /// structure, in addition to longitude, latitude and time.
  ///
  /// This parameter affects the order in which observations are processed and thus the final
  /// estimates of gross error probabilities, since the probability updates made when checking
  /// individual observation pairs are not commutative.
  oops::Parameter<bool> sortByPressure{"sort_by_pressure", false, this};

  /// Maximum total number of buddies of any observation.
  ///
  /// \note In the context of this parameter as well as the \c max_num_buddies_from_single_band
  /// \c max_num_buddies_with_same_station_id parameters, the number of buddies of any observation O
  /// is understood as the number of buddy pairs (O, O') where O' != O. This definition facilitates
  /// the buddy check implementation, but is an underestimate of the true number of buddies, since
  /// it doesn't take into account pairs of the form (O', O).
  oops::Parameter<int> maxTotalNumBuddies{"max_total_num_buddies", 15, this};

  /// Maximum number of buddies of any observation belonging to a single zonal band.
  ///
  /// See the note next to maxTotalNumBuddies.
  oops::Parameter<int> maxNumBuddiesFromSingleBand{"max_num_buddies_from_single_band", 10, this};

  /// Maximum number of buddies of any observation sharing that observation's station ID.
  ///
  /// See the note next to maxTotalNumBuddies.
  oops::Parameter<int> maxNumBuddiesWithSameStationId{
    "max_num_buddies_with_same_station_id", 5, this};

  /// Set to true to identify pairs of buddy observations using an algorithm reproducing exactly
  /// the algorithm used in Met Office's OPS system, but potentially skipping some valid buddy
  /// pairs.
  oops::Parameter<bool> useLegacyBuddyCollector{"use_legacy_buddy_collector", false, this};

  /// @}
  /// \name Parameters controlling gross error probability updates
  /// @{

  /// Encoding of the function mapping the latitude (in degrees) to the horizontal correlation scale
  /// (in km).
  ///
  /// The function is taken to be a linear interpolation of a series of (latitude, scale) points.
  /// The latitudes and scales at these points should be specified as keys and values of a
  /// JSON-style map. Owing to a bug in the eckit YAML parser, the keys must be enclosed in quotes.
  /// For example,
  ///
  ///   horizontal_correlation_scale: { "-90": 200, "90": 100 }
  ///
  /// encodes a function varying linearly from 200 km at the south pole to 100 km at the north pole.
  oops::Parameter<std::map<float, float>> horizontalCorrelationScaleInterpolationPoints{
    "horizontal_correlation_scale", {{-90.0f, 100.f}, {90.0f, 100.f}}, this};

  /// Temporal correlation scale.
  oops::Parameter<util::Duration> temporalCorrelationScale{"temporal_correlation_scale",
                                                           util::Duration("PT6H"), this};

  /// Parameter used to "damp" gross error probability updates using method 1 described in section
  /// 3.8 of the OPS Scientific Documentation Paper 2 to make the buddy check
  /// better-behaved in data-dense areas. See the reference above for the full description.
  oops::Parameter<double> dampingFactor1{"damping_factor_1", 1.0, this};

  /// Parameter used to "damp" gross error probability updates using method 2 described in section
  /// 3.8 of the OPS Scientific Documentation Paper 2 to make the buddy check
  /// better-behaved in data-dense areas. See the reference above for the full description.
  oops::Parameter<double> dampingFactor2{"damping_factor_2", 1.0, this};

  /// Non-divergence constraint. Used only for vector variables.
  oops::Parameter<double> nonDivergenceConstraint{"non_divergence_constraint", 1.0, this};

  /// @}
  /// \name Miscellaneous parameters
  /// @{

  /// Observations will be rejected if the gross error probability lies at or above this threshold.
  oops::Parameter<float> rejectionThreshold{"rejection_threshold", 0.5, this};

  /// Tracing information will be output for observations lying within any of the specified boxes.
  oops::Parameter<std::vector<LatLonBoxParameters>> tracedBoxes{"traced_boxes", {}, this};

  /// @}
};

}  // namespace ufo

#endif  // UFO_FILTERS_METOFFICEBUDDYCHECKPARAMETERS_H_
