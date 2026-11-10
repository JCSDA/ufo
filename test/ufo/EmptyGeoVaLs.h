/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UFO_EMPTYGEOVALS_H_
#define TEST_UFO_EMPTYGEOVALS_H_

#include <memory>
#include <utility>
#include <vector>

#include "ioda/ObsSpace.h"
#include "oops/base/Variables.h"
#include "oops/interface/SampledLocations.h"
#include "oops/util/DateTime.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsTraits.h"
#include "ufo/SampledLocations.h"

namespace ufo {
namespace test {

/// \brief Create a GeoVaLs object holding no variables and no locations, distributed like
/// \p obsspace.
///
/// Useful for tests of filters and predictors that require no GeoVaLs. Some filters treat
/// GeoVaLs with no locations as "no GeoVaLs supplied" (see e.g. ProfileDataHandler).
inline std::unique_ptr<GeoVaLs> makeEmptyGeoVaLs(const ioda::ObsSpace & obsspace) {
  const std::vector<float> lons, lats;
  const std::vector<util::DateTime> times;
  auto locs = std::make_unique<SampledLocations>(lons, lats, times, obsspace.distribution());
  return std::make_unique<GeoVaLs>(oops::SampledLocations<ObsTraits>(std::move(locs)),
                                   oops::Variables(), std::vector<size_t>());
}

}  // namespace test
}  // namespace ufo

#endif  // TEST_UFO_EMPTYGEOVALS_H_
