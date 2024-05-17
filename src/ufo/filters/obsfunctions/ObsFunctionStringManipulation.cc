/*
 * (C) Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsFunctionStringManipulation.h"

#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"

#include "oops/util/Logger.h"

namespace ufo {

constexpr char StringManipOptionParameterTraitsHelper::enumTypeName[];
constexpr util::NamedEnumerator<StringManipOption>
            StringManipOptionParameterTraitsHelper::namedValues[];
static ObsFunctionMaker<StringManipulation> makerStringManip_("StringManipulation");

// -----------------------------------------------------------------------------

StringManipulation::StringManipulation(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Check options
  options_.validateAndDeserialize(conf);

  // Create variable and add to invars_
  for (const Variable & var : options_.variable.value()) {
    invars_ += var;
  }
}

// -----------------------------------------------------------------------------

void StringManipulation::compute(const ObsFilterData & in,
                                 ioda::ObsDataVector<std::string> & out) const {
  // dimension
  const size_t nlocs = in.nlocs();

  // variables
  int startIndex;
  int cutLength;

  const std::string missing = util::missingValue<std::string>();
  ioda::ObsDataVector<std::string> varin(in.obsspace(), invars_.toOopsObsVariables());
  in.get(invars_[0], varin);

  // apply string operation to all locations
  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    switch (options_.stringManipOption.value()) {
      case StringManipOption::STRINGCUT: {
        if (options_.startIndex.value() != boost::none &&
            options_.cutLength.value() != boost::none) {
          // index to go from
          startIndex = options_.startIndex.value().get();
          // index to go to
          cutLength = options_.cutLength.value().get();
          // cut the input string from startIndex to a length cutLength
          out[0][iloc] = varin[0][iloc].substr(startIndex, cutLength);
        } else {
          throw eckit::Exception("startIndex and cutLength need to be set for stringcut");
        }
      }
    }
  }
}
// -----------------------------------------------------------------------------

const ufo::Variables & StringManipulation::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
