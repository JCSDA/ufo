/*
 * (C) Crown copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <sstream>

#include "ufo/utils/Constants.h"
#include "ufo/variabletransforms/Cal_Logarithm.h"

namespace ufo {
/******************************************************************************/
//  Cal_Logarithm
/******************************************************************************/
static TransformMaker<Cal_Logarithm> makerCal_Logarithm_("Logarithm");

Cal_Logarithm::Cal_Logarithm(
    const Parameters_ &options, const ObsFilterData &data,
    const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
    const std::shared_ptr<ioda::ObsDataVector<float>> &obserr)
    : TransformBase(options, data, flags, obserr),
      variable_(options.variable),
      group_(options.group) {
  if (options.base.value() != boost::none) {
    // If options.base.value() is boost::none, then base_ will be 0.0, which
    // indicates the natural logarithm.
    base_ = options.base.value().value_or(missingValueFloat);
    assert(base_ != missingValueFloat);  // This should never happen
    if (base_ <= 0.0 || base_ == 1.0) {
      std::ostringstream msg;
      msg << "Invalid logarithm base: " << base_;
      throw eckit::BadValue(msg.str(), Here());
    }
  }
  if (options.outputVariable.value() != boost::none) {
    outputVariable_ = options.outputVariable.value().value();
  }
  if (options.outputGroup.value() != boost::none) {
    outputGroup_ = options.outputGroup.value().value();
  }
}

/******************************************************************************/

void Cal_Logarithm::runTransform(const std::vector<bool> &apply) {
  oops::Log::trace() << "Retrieve variable" << std::endl;

  std::vector<float> variable;
  const size_t nlocs = obsdb_.nlocs();

  getObservation(group_, variable_, variable);

  if (variable.empty()) {
    std::ostringstream msg;
    msg << "Variable " << variable_ << " is empty";
    throw eckit::BadValue(msg.str(), Here());
  }

  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    // Set to missing if the data have been excluded by the where statement,
    // are missing or calculation of the logarithm is not possible.
    if (!apply[iloc] || variable[iloc] == missingValueFloat ||
        variable[iloc] <= 0.0) {
      variable[iloc] = missingValueFloat;
    } else if (base_ == 0.0) {  // Natural logarithm (base e)
      variable[iloc] = std::log(variable[iloc]);
    } else if (base_ == 10.0) {
      variable[iloc] = std::log10(variable[iloc]);
    } else if (base_ == 2.0) {
      variable[iloc] = std::log2(variable[iloc]);
    } else {
      variable[iloc] = std::log(variable[iloc]) / std::log(base_);
    }
  }

  // Where the group is not specified, the derived version of the input group
  // name is used.
  if (outputGroup_.empty() && !outputVariable_.empty()) {
    obsdb_.put_db(getDerivedGroup(group_), outputVariable_, variable);
  } else if (!outputGroup_.empty() && outputVariable_.empty()) {
    obsdb_.put_db(outputGroup_, variable_, variable);
  } else if (!outputGroup_.empty() && !outputVariable_.empty()) {
    obsdb_.put_db(outputGroup_, outputVariable_, variable);
  } else {
    obsdb_.put_db(getDerivedGroup(group_), variable_, variable);
  }
}
}  // namespace ufo
