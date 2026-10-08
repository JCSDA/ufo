/*
 * (C) Copyright 2019 UK Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */


#ifndef UFO_OPERATORS_TIMEOPER_OBSTIMEOPERUTIL_H_
#define UFO_OPERATORS_TIMEOPER_OBSTIMEOPERUTIL_H_

#include <algorithm>
#include <ostream>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Locations.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace ufo {

class GeoVaLs;
class ObsTimeOperParameters;
struct ObsTraits;

std::vector<std::vector<float>> timeWeightCreate(const ioda::ObsSpace & odb_,
                                                 const ObsTimeOperParameters & parameters);

/// Two interpolation paths per location, at the model states before and after the observation.
oops::Locations<ObsTraits> timeOperLocations(const ioda::ObsSpace & odb,
                                             const util::Duration & windowSub);

/// Store in the reduced format the time-weighted sum of the two profiles sampling each location.
void timeInterpolate(GeoVaLs & gv, const std::vector<std::vector<float>> & timeWeights);

/// Adjoint of timeInterpolate: add the weighted reduced profiles of `gvad` to the sampled `gv`.
void timeInterpolateAD(GeoVaLs & gv, const GeoVaLs & gvad,
                       const std::vector<std::vector<float>> & timeWeights);

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_TIMEOPER_OBSTIMEOPERUTIL_H_
