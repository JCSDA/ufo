/*
 * (C) British Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_PROFILEUNFLAGOBSCHECK_H_
#define UFO_FILTERS_PROFILEUNFLAGOBSCHECK_H_

#include <map>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/NumericConstraints.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

/// Parameters controlling the operation of the ProfileUnFlagObsCheck filter.
class ProfileUnFlagObsCheckParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(ProfileUnFlagObsCheckParameters, FilterParametersBase)

 public:
  /// The filter will unflag observations where surrounding values match within
  /// a tolerance
  oops::RequiredParameter<float> aTol{"absolute tolerance", this};
  /// The function is taken to be a linear interpolation of a series of
  /// (vertical coordinate, scale) points. The position and scales at these
  /// points should be specified as keys and values of a JSON-style map. Owing
  /// to a bug in the eckit YAML parser, the keys must be enclosed in quotes.
  /// For example,
  ///
  /// \code{.yaml}
  ///   vertical tolerance scale: { "0": 1, "2": 0.5 }
  /// \endcode
  ///
  /// encodes a function varying linearly vertically.
  oops::OptionalParameter<std::map<float, float>> verticalToleranceScaleInterpolationPoints{
      "vertical tolerance scale", this};
  /// Name of vertical coordinate
  oops::Parameter<Variable> verticalCoord{"vertical coordinate",
                                          "Name of vertical coordinate variable",
                                          Variable{"ObsValue/depthBelowWaterSurface"},
                                          this};
};

/// \brief ProfileUnFlagObsCheck: Unflag observations in a profile that seem correct
///
/// Uses the record number in obsdb_ to identify which observations belong to a
/// given profile (all members of a profile must share the same record number).
/// Each observation in a profile is compared to those above and below. If both
/// of these are unflagged and match the observation to within a tolerance,
/// then the observation is unflagged. If the observation is the first or the
/// last in the profile than a match with only the single adjacent observation
/// is sufficient for unflagging.
///
/// The tolerance is set using an absolute value and an optional piecewise
/// interpolation to scale the tolerance based on the vertical coordinate.
///
/// Example:
///
/// \code{.yaml}
///  - filter: Profile Unflag Observations Check
///    filter variables:
///    - name: ObsValue/waterTemperature
///    absolute tolerance:  1.0
///    vertical tolerance scale: { "0": 1, "2": 0.5 }
///    vertical coordinate: ObsValue/depthBelowWaterSurface
/// \endcode
///

class ProfileUnFlagObsCheck : public FilterBase,
                           private util::ObjectCounter<ProfileUnFlagObsCheck> {
 public:
  /// The type of parameters accepted by the constructor of this filter.
  /// This typedef is used by the FilterFactory.
  typedef ProfileUnFlagObsCheckParameters Parameters_;

  static const std::string classname() {return "ufo::ProfileUnFlagObsCheck";}

  ProfileUnFlagObsCheck(ioda::ObsSpace &, const Parameters_ &,
                     std::shared_ptr<ioda::ObsDataVector<int> >,
                     std::shared_ptr<ioda::ObsDataVector<float> >);
  ~ProfileUnFlagObsCheck();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::pass;}
  Parameters_ parameters_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_PROFILEUNFLAGOBSCHECK_H_
