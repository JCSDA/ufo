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

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/ObsOperatorBase.h"

#include "ufo/profile/ObsProfileAverageData.h"
#include "ufo/profile/ObsProfileAverageParameters.h"

/// Forward declarations
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
/// intersect each model pressure level. The intersections are found by stepping through the
/// observation locations from the lowest-altitude value upwards. For each model level,
/// the location of the observation whose pressure is larger than, and closest to, the model
/// pressure is recorded. The \p<vertical coordinate> parameter controls the model pressure
/// GeoVaLs that are used in this procedure.
/// If there are no observations in a model level, which can occur for (e.g.) sondes reporting in
/// low-frequency TAC format, the location corresponding to the last filled level is used.
/// (If there are some model levels closer to the surface than the lowest-altitude observation,
/// the location of the lowest observation is used for these levels.)
///
/// This procedure is iterated multiple times in order to account for the fact that model pressures
/// can be slanted close to the Earth's surface.
/// The number of iterations is configured with the \p<number of intersection iterations> parameter.
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
///      allocate companion records with length: 71
///
/// (where 71 can be replaced by the length of the air_pressure_levels GeoVaL).
/// The H(x) values are placed in the extended section of the ObsSpace.
/// Note that, unlike what may be expected for an observation operator, averaging of the model
/// values across each layer is not performed; a single model value is used in each case.
/// This follows what is used in OPS. Alternatives could be considered in the future.
///
/// A comparison with OPS is be performed if the option \p<compare with OPS> is set to true.
/// This checks values of the locations and pressure values associated with the slant path.
/// All other comparisons are performed with the standard 'vector ref' option in the yaml file.
///
/// This operator also accepts an optional `variables` parameter, which controls which ObsSpace
/// variables will be simulated. This option should only be set if this operator is used as a
/// component of the Composite operator. If `variables` is not set, the operator will simulate
/// all ObsSpace variables. Please see the documentation of the Composite operator for further
/// details.
class ObsProfileAverage : public ObsOperatorBase,
  private util::ObjectCounter<ObsProfileAverage> {
 public:
  static const std::string classname() {return "ufo::ObsProfileAverage";}
  typedef ioda::ObsDataVector<int> QCFlags_t;
  typedef ObsProfileAverageParameters Parameters_;

  ObsProfileAverage(const ioda::ObsSpace &, const Parameters_ &);
  ~ObsProfileAverage() override;

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t &) const override;

  const oops::Variables & requiredVars() const override { return requiredVars_; }

  oops::ObsVariables simulatedVars() const override { return data_.simulatedVars(); }

 private:
  void print(std::ostream &) const override;

 private:
  /// ObsSpace.
  const ioda::ObsSpace& odb_;

  /// Data handler for the ProfileAverage operator and TL/AD code.
  ObsProfileAverageData data_;

  /// Required variables.
  oops::Variables requiredVars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_PROFILE_OBSPROFILEAVERAGE_H_
