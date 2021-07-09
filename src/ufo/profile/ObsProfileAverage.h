/*
 * (C) Copyright 2021 UK Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_OBSPROFILEAVERAGE_H_
#define UFO_PROFILE_OBSPROFILEAVERAGE_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/ObsOperatorBase.h"

#include "ufo/profile/ObsProfileAverageData.h"
#include "ufo/profile/ObsProfileAverageParameters.h"

/// Forward declarations
namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class Locations;
  class ObsDiagnostics;

/// \brief Observation operator for profile averaging.
///
/// This observation operator produces H(x) vectors which correspond to vertically-averaged
/// profiles. The algorithm determines the locations at which reported-level profiles
/// intersect each model pressure level (based on the air_pressure_levels GeoVaL).
/// This is done by stepping through the observation locations from the lowest-altitude value
/// upwards. For each model level, the location of the observation whose pressure is larger than,
/// and closest to, the model pressure is recorded.
///
/// If there are no observations in a model level, which can occur for (e.g.) sondes reporting in
/// low-frequency TAC format, the location corresponding to the last filled level is used.
/// (If there are some model levels closer to the surface than the lowest-altitude observation,
/// the location of the lowest observation is used for these levels.)
///
/// This procedure is iterated multiple times in order to account for the fact that model pressures
/// can be slanted close to the Earth's surface.
/// The number of iterations is configured with the \p numIntersectionIterations parameter.
///
/// Having obtained the profile boundaries, values of model pressure and any simulated variables
/// are obtained as in the locations that were determined in the procedure above.
/// This produces a single column of model values which are used as the H(x) variable.
/// In essence, this operator converts a set of GeoVaLs to what is referred to as a 'CX column'
/// in OPS terminology.
///
/// In order for this operator to work correctly the ObsSpace must have been extended as in
/// the following yaml snippet:
///
/// - obs space:
///    extension:
///      average profiles onto model levels: 71
///
/// (where 71 can be replaced by the length of the air_pressure_levels GeoVaL).
/// The H(x) values are placed in the extended section of the ObsSpace.
/// Note that, unlike what may be expected for an observation operator, averaging of the model
/// values across each layer is not performed; a single model value is used in each case.
/// This follows what is used in OPS. Alternatives could be considered in the future.
///
/// A comparison with OPS is be performed if the option \p compareWithOPS is set to true.
/// This checks values of the locations and pressure values associated with the slant path.
/// All other comparisons are performed with the standard 'vector ref' option in the yaml file.
class ObsProfileAverage : public ObsOperatorBase,
  private util::ObjectCounter<ObsProfileAverage> {
 public:
  static const std::string classname() {return "ufo::ObsProfileAverage";}

  ObsProfileAverage(const ioda::ObsSpace &, const eckit::Configuration &);
  ~ObsProfileAverage() override;

  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &) const override;

  const oops::Variables & requiredVars() const override { return data_.requiredVars(); }

 private:
  void print(std::ostream &) const override;

 private:
  /// ObsSpace.
  const ioda::ObsSpace& odb_;

  /// Data handler for the ProfileAverage operator and TL/AD code.
  ObsProfileAverageData data_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_PROFILE_OBSPROFILEAVERAGE_H_
