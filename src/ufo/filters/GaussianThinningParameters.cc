/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/GaussianThinningParameters.h"

#include <utility>

namespace ufo {
constexpr char DistanceNormParameterTraitsHelper::enumTypeName[];
constexpr util::NamedEnumerator<DistanceNorm> DistanceNormParameterTraitsHelper::namedValues[];

void GaussianThinningParameters::deserialize(util::CompositePath &path,
                                             const eckit::Configuration &config) {
  oops::Parameters::deserialize(path, config);

  if (opsCompatibilityMode) {
    if (roundHorizontalBinCountToNearest.value() != boost::none &&
        *roundHorizontalBinCountToNearest.value() == false)
      throw eckit::UserError(
            path.path() + ": round_horizontal_bin_count_to_nearest must not be set to false when "
                          "ops_compatibility_mode is set to true", Here());
    if (partitionLongitudeBinsUsingMesh.value() != boost::none &&
        *partitionLongitudeBinsUsingMesh.value() == false)
      throw eckit::UserError(
            path.path() + ": partition_longitude_bins_using_mesh must not be set to false when "
                          "ops_compatibility_mode is set to true", Here());
    if (defineMeridian20000km.value() != boost::none &&
        *defineMeridian20000km.value() == false)
      throw eckit::UserError(
            path.path() + ": define_meridian_2000_km must not be set to false when "
                          "ops_compatibility_mode is set to true", Here());
    if (distanceNorm.value() != boost::none &&
        *distanceNorm.value() == DistanceNorm::GEODESIC)
      throw eckit::UserError(
            path.path() + ": distance_norm must not be set to 'geodesic' when "
                          "ops_compatibility_mode is set to true", Here());
  }
}

}  // namespace ufo
