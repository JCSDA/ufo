/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/filters/ObsBoundsCheck.h"

#include <algorithm>
#include <set>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/QCflags.h"
#include "ufo/utils/PrimitiveVariables.h"
#include "ufo/utils/StringUtils.h"

namespace ufo {

namespace {

// -----------------------------------------------------------------------------

/// Set each element of \p flagged to true if the corresponding element of \p apply is true
/// and the corresponding element of \p testValues is a non-missing value lying outside the
/// closed interval [\p minValue, \p maxValue].
void flagWhereOutOfBounds(const std::vector<bool> & apply,
                          const std::vector<float> & testValues,
                          const float minValue,
                          const float maxValue,
                          const bool treatMissingAsOutOfBounds,
                          std::vector<bool> &flagged) {
  const size_t nlocs = testValues.size();
  ASSERT(apply.size() == nlocs);
  ASSERT(flagged.size() == nlocs);

  const float missing = util::missingValue<float>();
  for (size_t i = 0; i < nlocs; ++i) {
    if (apply[i]) {
      if (testValues[i] != missing) {
        if (minValue != missing && testValues[i] < minValue)
          flagged[i] = true;
        if (maxValue != missing && testValues[i] > maxValue)
          flagged[i] = true;
      } else {
        if (treatMissingAsOutOfBounds)
          flagged[i] = true;
      }
    }
  }
}

}  // namespace

// -----------------------------------------------------------------------------

ObsBoundsCheck::ObsBoundsCheck(ioda::ObsSpace & obsdb, const Parameters_ & parameters,
                               std::shared_ptr<ioda::ObsDataVector<int> > flags,
                               std::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, parameters, flags, obserr), parameters_(parameters)
{
  if (parameters_.testVariables.value() != boost::none) {
    for (const Variable & var : *parameters_.testVariables.value())
      allvars_ += var;
  }
  oops::Log::debug() << "ObsBoundsCheck: config (constructor) = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

ObsBoundsCheck::~ObsBoundsCheck() {}

// -----------------------------------------------------------------------------

void ObsBoundsCheck::applyFilter(const std::vector<bool> & apply,
                                 const Variables & filtervars,
                                 std::vector<std::vector<bool>> & flagged) const {
  // Find the variables that should be tested. Use the variables specified in the 'test variables'
  // option if present, otherwise the filter variables.
  ufo::Variables testvars;
  if (parameters_.testVariables.value() != boost::none) {
    for (const Variable & var : *parameters_.testVariables.value())
      testvars += var;
  } else {
    testvars += ufo::Variables(filtervars, "ObsValue");
  }
  if (!testvars)
    throw eckit::UserError("ObsBoundsCheck: The list of test variables is empty", Here());

  oops::Log::debug() << "ObsBoundsCheck: filtering " << filtervars << " with "
                     << testvars << std::endl;

  // Retrieve the bounds.
  const float missing = util::missingValue<float>();
  const float vmin = parameters_.minvalue.value().value_or(missing);
  const float vmax = parameters_.maxvalue.value().value_or(missing);

  // Determine the mode of operation.
  const bool flagAllFilterVarsIfAnyTestVarOutOfBounds =
      parameters_.testVariables.value() != boost::none &&
      (parameters_.flagAllFilterVarsIfAnyTestVarOutOfBounds.value() ||
       testvars.nvars() == 1);
  const bool onlyTestGoodFilterVarsForFlagAllFilterVars =
      parameters_.onlyTestGoodFilterVarsForFlagAllFilterVars;
  const bool treatMissingAsOutOfBounds = parameters_.treatMissingAsOutOfBounds;

  // Do the actual work.
  if (flagAllFilterVarsIfAnyTestVarOutOfBounds) {
    // If using test only good filter vars then the number of filter
    // vars must equal the number of test vars
    if (onlyTestGoodFilterVarsForFlagAllFilterVars && filtervars.nvars() != testvars.nvars())
      throw eckit::UserError("The number of 'primitive' (single-channel) test variables must "
                             "match that of 'primitive' filter variables when using the "
                             "'test only filter variables with passed qc when flagging all "
                             "filter variables' option.");
    ASSERT(filtervars.nvars() == flagged.size());
    // Convert to an oops variables which contains a simple list of variables
    // expanding out the channels and therefore it matches the testvars list
    const oops::ObsVariables filtervarslist = filtervars.toOopsObsVariables();
    // Loop over all channels of all test variables and record all locations where any of these
    // channels is out of bounds.
    std::vector<bool> anyTestVarOutOfBounds(obsdb_.nlocs(), false);
    size_t ifiltervar = 0;
    for (PrimitiveVariable singleChannelTestVar : PrimitiveVariables(testvars, data_)) {
      std::vector<bool> testAtLocations = apply;
      if (onlyTestGoodFilterVarsForFlagAllFilterVars) {
        for (size_t iloc=0; iloc < testAtLocations.size(); iloc++)
          if ((*flags_)[filtervarslist[ifiltervar]][iloc] != QCflags::pass) {
            testAtLocations[iloc] = false;
          }
      }
      const std::vector<float> & testValues = singleChannelTestVar.values();
      flagWhereOutOfBounds(testAtLocations, testValues, vmin, vmax, treatMissingAsOutOfBounds,
                           anyTestVarOutOfBounds);
      ifiltervar++;
    }
    // Copy these flags to the flags of all filtered variables.
    for (std::vector<bool> &f : flagged)
      f = anyTestVarOutOfBounds;
  } else {
    if (filtervars.nvars() != testvars.nvars())
      throw eckit::UserError("The number of 'primitive' (single-channel) test variables must match "
                             "that of 'primitive' filter variables unless the 'flag all filter "
                             "variables if any test variable is out of bounds' option is set");
    // Loop over all channels of all test variables and for each locations where that channel is out
    // of bounds, flag the corresponding filter variable channel.
    ASSERT(filtervars.nvars() == flagged.size());
    size_t ifiltervar = 0;
    for (PrimitiveVariable singleChannelTestVar : PrimitiveVariables(testvars, data_)) {
      const std::vector<float> & testValues = singleChannelTestVar.values();
      flagWhereOutOfBounds(apply, testValues, vmin, vmax, treatMissingAsOutOfBounds,
                           flagged[ifiltervar++]);
    }
  }
}

// -----------------------------------------------------------------------------

void ObsBoundsCheck::print(std::ostream & os) const {
  os << "ObsBoundsCheck: config = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
