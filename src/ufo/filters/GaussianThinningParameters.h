/*
 * (C) Copyright 2019 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_GAUSSIANTHINNINGPARAMETERS_H_
#define UFO_FILTERS_GAUSSIANTHINNINGPARAMETERS_H_

#include <string>

#include "eckit/exception/Exceptions.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "ufo/utils/Constants.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace eckit {
  class Configuration;
}

namespace ufo {

enum class DistanceNorm {
  GEODESIC, MAXIMUM
};

}  // namespace ufo

namespace oops {

template <>
struct ParameterTraits<ufo::DistanceNorm> {
  static boost::optional<ufo::DistanceNorm> get(const eckit::Configuration &config,
                                                const std::string& name) {
    std::string value;
    if (config.get(name, value)) {
      if (value == "geodesic")
        return ufo::DistanceNorm::GEODESIC;
      if (value == "maximum")
        return ufo::DistanceNorm::MAXIMUM;
      throw eckit::BadParameter("Bad conversion from std::string '" + value + "' to DistanceNorm",
                                Here());
    } else {
      return boost::none;
    }
  }
};

}  // namespace oops

namespace ufo {

/// \brief Options controlling the operation of the Gaussian_Thinning filter.
class GaussianThinningParameters : public oops::Parameters {
 public:
  // Horizontal grid

  /// Cell size (in km) along the meridians. Thinning in the horizontal direction is disabled if
  /// this parameter is negative.
  // FIXME(wsmigaj): for consistency with vertical thinning, I think it would be better to
  // interpret absence of the horizontal_mesh setting as "don't do horizontal thinning"
  // rather than using an arbitrary default grid size. Leaving things as they are for backward
  // compatibility, for now.
  oops::Parameter<float> horizontalMesh{"horizontal_mesh", defaultHorizontalMesh(), this};
  /// True to use a reduced grid, with fewer cells at high latitudes.
  /// False to use a regular grid, with the same number of cells at all latitudes.
  oops::Parameter<bool> useReducedHorizontalGrid{"use_reduced_horizontal_grid", true, this};
  /// True to set the number of zonal bands so that the band width is as close as possible to
  /// \c horizontal_mesh, and the number of bins in each zonal band so that the bin width in the
  /// zonal direction is as close as possible to that in the meridional direction.
  /// False to set the number of zonal bands so that the band width is as small as possible, but
  /// no smaller than \c horizontal_mesh, and the bin width in the zonal direction is as small as
  /// possible, but no smaller than in the meridional direction.
  oops::Parameter<bool> roundHorizontalBinCountToNearest{
    "round_horizontal_bin_count_to_nearest", false, this};

  // Vertical grid

  /// Cell size (in Pa) in the vertical direction. Thinning in the vertical direction is disabled
  /// if this parameter is not specified or negative.
  oops::Parameter<float> verticalMesh{"vertical_mesh", -1.0f, this};
  /// Lower bound of the pressure interval split into cells of size \c vertical_mesh.
  oops::Parameter<float> verticalMin{"vertical_min", 100.0f, this};
  /// Upper bound of the pressure interval split into cells of size \c vertical_mesh.
  /// This parameter is rounded upwards to the nearest multiple of \c vertical_mesh starting from
  /// \c vertical_min.
  oops::Parameter<float> verticalMax{"vertical_max", 110000.0f, this};

  // Temporal grid

  /// Cell size in the temporal direction. Temporal thinning is disabled if this this parameter is
  /// not specified or set to 0.
  oops::OptionalParameter<util::Duration> timeMesh{"time_mesh", this};
  /// Lower bound of the time interval split into cells of size \c time_mesh. Temporal thinning is
  /// disabled if this parameter is not specified.
  oops::OptionalParameter<util::DateTime> timeMin{"time_min", this};
  /// Upper bound of the time interval split into cells of size \c time_mesh.
  /// This parameter is rounded upwards to the nearest multiple of \c time_mesh starting from
  /// \c time_min. Temporal thinning is disabled if this parameter is not specified.
  oops::OptionalParameter<util::DateTime> timeMax{"time_max", this};

  // Observation categories

  /// Variable storing integer-valued IDs associated with observations. Observations belonging
  /// to different categories are thinned separately.
  oops::OptionalParameter<Variable> categoryVariable{"category_variable", this};

  // Selection of observations to retain

  /// Variable storing observation priorities. Among all observations in a cell, only those with
  /// the highest priority are considered as candidates for retaining. If not specified, all
  /// observations are assumed to have equal priority.
  oops::OptionalParameter<Variable> priorityVariable{"priority_variable", this};

  /// Determines which of the highest-priority observations lying in a cell is retained.
  /// Allowed values:
  /// - \c geodesic: retain the observation closest to the cell centre in the horizontal direction
  ///   (air pressure and time are ignored)
  /// - \c maximum: retain the observation lying furthest from the cell's bounding box in the
  ///   system of coordinates in which the cell is a unit cube (all dimensions along which thinning
  ///   is enabled are taken into account).
  oops::Parameter<DistanceNorm> distanceNorm{"distance_norm", DistanceNorm::GEODESIC, this};

 private:
  static float defaultHorizontalMesh() {
    return static_cast<float>(2 * M_PI * Constants::mean_earth_rad / 360.0);
  }
};

}  // namespace ufo

#endif  // UFO_FILTERS_GAUSSIANTHINNINGPARAMETERS_H_
