/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_DATETIMEOFFSET_H_
#define UFO_FILTERS_OBSFUNCTIONS_DATETIMEOFFSET_H_

#include <string>
#include <vector>

#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"

namespace ufo {

enum class OffsetUnit {
  SECONDS, MINUTES, HOURS
};

struct OffsetUnitParameterTraitsHelper {
  typedef OffsetUnit EnumType;
  static constexpr char enumTypeName[] = "OffsetUnit";
  static constexpr util::NamedEnumerator<OffsetUnit> namedValues[] = {
    { OffsetUnit::SECONDS, "seconds" },
    { OffsetUnit::MINUTES, "minutes" },
    { OffsetUnit::HOURS, "hours" }
  };
};

}  // namespace ufo

namespace oops {
template <>
struct ParameterTraits<ufo::OffsetUnit> :
  public EnumParameterTraits<ufo::OffsetUnitParameterTraitsHelper>
{};
}

namespace ufo {

/// \brief Options controlling the DateTimeOffset ObsFunction
class DateTimeOffsetParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(DateTimeOffsetParameters, Parameters)

 public:
  oops::RequiredParameter<std::string> offset_name
    {"offset variable name",
     "Name of the offset variable.",
     this};

  oops::RequiredParameter<OffsetUnit> offset_unit
    {"offset unit",
     "Name of the offset unit. Valid options: seconds, minutes, hours.",
     this};
};

// -----------------------------------------------------------------------------

/// \brief Add time offsets (in seconds, minutes or hours) to dateTime@MetaData.
///
/// This function is used to add a location-dependent offset to dateTime@MetaData.
/// If the offset at a particular location is missing then
/// the corresponding datetime is not modified.
class DateTimeOffset : public ObsFunctionBase<util::DateTime> {
 public:
  explicit DateTimeOffset(const eckit::LocalConfiguration &);
  ~DateTimeOffset();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<util::DateTime> &) const;
  const ufo::Variables & requiredVariables() const;

 private:
  // Add offset to DateTime at each location (unless offset is missing).
  template <typename T>
  void applyOffsets(const ObsFilterData & in,
                    ioda::ObsDataVector<util::DateTime> & out) const {
    const T missing = util::missingValue(missing);
    const size_t nlocs = in.obsspace().nlocs();

    // Get datetime.
    std::vector <util::DateTime> datetimes(nlocs);
    in.get(Variable("MetaData/dateTime"), datetimes);

    // Get offset variable.
    std::vector <T> offsets(nlocs);
    in.get(Variable(options_.offset_name.value()), offsets);

    // Get offset multiplier, which depends on the offset unit.
    int offsetmult = 1;
    if (options_.offset_unit == OffsetUnit::MINUTES) offsetmult = 60;
    else if (options_.offset_unit == OffsetUnit::HOURS) offsetmult = 3600;

    // Loop through locations and apply any non-missing offsets to datetime.
    for (size_t jloc = 0; jloc < nlocs; ++jloc) {
      const T offset = offsets[jloc];
      // If the offset is missing do not modify the datetime.
      if (offset == missing)
        out[0][jloc] = datetimes[jloc];
      else
        out[0][jloc] = datetimes[jloc] +
          util::Duration(static_cast<int64_t>(offset * offsetmult));
    }
  }

  DateTimeOffsetParameters options_;
  ufo::Variables invars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_DATETIMEOFFSET_H_
