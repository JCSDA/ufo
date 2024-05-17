/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "eckit/utils/StringTools.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/utils/PrimitiveVariables.h"

namespace ufo {

void PrimitiveVariablesIterator::loadCurrentVariable() {
  if (variableIndex_ < variables_.size()) {
    const Variable &variable = variables_[variableIndex_];
    vector_.reset(new ioda::ObsDataVector<float>(data_.obsspace(), variable.toOopsObsVariables()));
    if (eckit::StringTools::endsWith(variable.group(), "ObsFunction")) {
      data_.get(variable, *vector_);
    } else {
      for (size_t i = 0; i < variable.size(); ++i) {
        data_.get(variable[i], (*vector_)[i]);
      }
    }
  }
}

}  // namespace ufo
