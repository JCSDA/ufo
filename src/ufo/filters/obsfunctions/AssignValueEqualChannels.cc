/* -----------------------------------------------------------------------------
 * (C) British Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * -----------------------------------------------------------------------------
 */

/* -----------------------------------------------------------------------------
 * Function to assign a value if the input matches a test value. If not, then
 * assign an alternative value.
 * -----------------------------------------------------------------------------
 */
#include "ufo/filters/obsfunctions/AssignValueEqualChannels.h"

#include <math.h>
#include <algorithm>
#include <set>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variable.h"

#include "eckit/exception/Exceptions.h"

namespace ufo {

static ObsFunctionMaker<AssignValueEqualChannels>
  makerAssignValueEqualChannels_("AssignValueEqualChannels");

/* -----------------------------------------------------------------------------
 * Constructor: Get the options and specify that the input variable must be read
 * -----------------------------------------------------------------------------
 */
AssignValueEqualChannels::AssignValueEqualChannels(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Get options into a local variable
  options_.deserialize(conf);

  // Create variable and add to invars_
  invars_ += options_.variable.value();
}

// -----------------------------------------------------------------------------

AssignValueEqualChannels::~AssignValueEqualChannels() {}

/* -----------------------------------------------------------------------------
 * Perform the computation.
 * Check that the channels given matches the number of variables to output.
 * Look through input, and assign output values if the input matches the test
 * value.
 * -----------------------------------------------------------------------------
 */
void AssignValueEqualChannels::compute(const ObsFilterData & in,
                        ioda::ObsDataVector<float> & out) const {
  const size_t nlocs = in.nlocs();
  const size_t nchans = out.nvars();
  const std::vector<int> channels = options_.variable.value().channels();

  if (nchans == 1 && channels.size() == 0) {
    // All good, do nothing
  } else if (nchans != channels.size()) {
    std::string errString = "AssignValueEqualChannels: mismatch between channels (" +
      std::to_string(channels.size()) + ") and number of variables (" +
      std::to_string(nchans) + ")";
    throw eckit::BadValue(errString);
  }

  size_t ivar = 0;
  for (size_t ichan=0; ichan < nchans ; ichan++) {
    std::vector<int> qcflags;
    in.get(Variable(options_.variable.value().fullName(), channels)[ichan], qcflags);
    for (size_t jj = 0; jj < nlocs; ++jj) {
      if (qcflags[jj] == options_.testValue.value()) {
          out[ivar][jj] = options_.assignEqual.value();
      } else {
          out[ivar][jj] = options_.assignNotEqual.value();
      }
    }
    ivar++;
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & AssignValueEqualChannels::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
