/*
 * (C) Crown copyright 2026 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_LOCATIONTOOBSSPACE_H_
#define UFO_FILTERS_OBSFUNCTIONS_LOCATIONTOOBSSPACE_H_

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

enum class LocationToObsSpaceLabelType {
  GlobalLocation, InRecordLocation
};

struct LocationToObsSpaceLabelTypeParameterTraitsHelper {
  typedef LocationToObsSpaceLabelType EnumType;
  static constexpr char enumTypeName[] = "LocationToObsSpaceLabelType";
  static constexpr util::NamedEnumerator<LocationToObsSpaceLabelType> namedValues[] = {
    { LocationToObsSpaceLabelType::GlobalLocation, "global location" },
    { LocationToObsSpaceLabelType::InRecordLocation, "in-record location" }
  };
};

}  // namespace ufo

namespace oops {

template <>
struct ParameterTraits<ufo::LocationToObsSpaceLabelType> :
    public EnumParameterTraits<ufo::LocationToObsSpaceLabelTypeParameterTraitsHelper>
{};

}  // namespace oops

namespace ufo {

/// \brief Parameters controlling the operation of the LocationToObsSpace ObsFunction.
class LocationToObsSpaceParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(LocationToObsSpaceParameters, Parameters)

 public:
  /// label type, either "in-record location" or "global location" (default).
  oops::Parameter<LocationToObsSpaceLabelType> label_type{"label type",
    LocationToObsSpaceLabelType::GlobalLocation, this};
};

/// \brief Write a location index to a variable which can be saved to the ObsSpace.
/// If `label type` is "global location" (default), write the global location index.
/// If `label type` is "in-record location", write the 0-based index within each
/// record (requires record grouping).
class LocationToObsSpace : public ObsFunctionBase<int> {
 public:
  explicit LocationToObsSpace(const eckit::LocalConfiguration &);
  ~LocationToObsSpace();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<int> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  LocationToObsSpaceParameters options_;
  ufo::Variables invars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_LOCATIONTOOBSSPACE_H_
