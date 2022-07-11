/*
 * (C) Copyright 2021 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/DateTimeOffset.h"

#include "ioda/ObsDataVector.h"

namespace ufo {

static ObsFunctionMaker<DateTimeOffset>
                       makerDateTimeOffset_("DateTimeOffset");

constexpr char OffsetUnitParameterTraitsHelper::enumTypeName[];
constexpr util::NamedEnumerator<OffsetUnit>
  OffsetUnitParameterTraitsHelper::namedValues[];

// -----------------------------------------------------------------------------

DateTimeOffset::DateTimeOffset(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Validate and deserialize options
  options_.validateAndDeserialize(conf);

  invars_ += Variable("MetaData/dateTime");
  invars_ += Variable(options_.offset_name.value());
}

// -----------------------------------------------------------------------------

DateTimeOffset::~DateTimeOffset() {}

// -----------------------------------------------------------------------------

void DateTimeOffset::compute(const ObsFilterData & in,
                             ioda::ObsDataVector<util::DateTime> & out) const {
  oops::Log::trace() << "DateTimeOffset::compute started" << std::endl;

  // Apply offsets to datetimes.
  // todo(ctgh): add Integer_64 when it becomes available.
  if (in.dtype(Variable(options_.offset_name.value())) == ioda::ObsDtype::Integer)
    applyOffsets<int>(in, out);
  else if (in.dtype(Variable(options_.offset_name.value())) == ioda::ObsDtype::Float)
    applyOffsets<float>(in, out);
  else
    throw eckit::BadParameter("Offset variable has incorrect type", Here());

  oops::Log::trace() << "DateTimeOffset::compute finished" << std::endl;
}

// -----------------------------------------------------------------------------

const ufo::Variables & DateTimeOffset::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
