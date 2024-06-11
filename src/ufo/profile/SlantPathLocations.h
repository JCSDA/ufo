/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_PROFILE_SLANTPATHLOCATIONS_H_
#define UFO_PROFILE_SLANTPATHLOCATIONS_H_

#include <string>
#include <vector>

namespace ioda {
  class ObsSpace;
}

namespace ufo {
  class GeoVaLs;
}

namespace ufo {
  /// Get slant path locations. This determines, for each model level, the location that
  /// corresponds to the intersection of the observed profile with that level.
  ///
  /// \param odb
  ///   ObsSpace.
  /// \param gv
  ///   GeoVaLs.
  /// \param locs
  ///   All locations in the profile.
  /// \param modelVerticalCoord
  ///   Name of the vertical coordinate used in the model.
  /// \param obsVerticalCoord
  ///   The full name (e.g. MetaData/air_pressure) of the observed vertical coordinate.
  /// \param itermax
  ///   Maximum number of interations that will be used to find the intersections
  ///   between observed pressures and model levels.
  ///
  /// \returns A vector of the slant path locations.
  std::vector<std::size_t> getSlantPathLocations(const ioda::ObsSpace & odb,
                                                 const GeoVaLs & gv,
                                                 const std::vector<std::size_t> & locs,
                                                 const std::string & obsVerticalCoord,
                                                 const oops::Variable & modelVerticalCoord,
                                                 const int itermax = 3);
}  // namespace ufo

#endif  // UFO_PROFILE_SLANTPATHLOCATIONS_H_
