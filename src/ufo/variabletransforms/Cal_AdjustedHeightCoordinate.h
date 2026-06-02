/*
 * (C) Crown copyright 2021, Met Office
 * (C) Copyright 2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "ufo/variabletransforms/TransformBase.h"

namespace ufo {

// -------------------------------------------------------------------------------------------------

  /*!
  * \brief Calculate the adjusted height coordinate with station elevation offset
  *
  * \details  Within the Variable Transforms filter, apply the transform "AdjustedHeightCoordinate"
  *  in order to calculate the the the new height coordinate taking station elevation into account.
  */

  class Cal_AdjustedHeightCoordinate : public TransformBase {
   public:
    Cal_AdjustedHeightCoordinate(const GenericVariableTransformParameters &options,
              const ObsFilterData &data,
              const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
              const std::shared_ptr<ioda::ObsDataVector<float>> &obserr);
    // Run variable conversion
    void runTransform(const std::vector<bool> &apply) override;
    Variables requiredVariables() const override { return gvals_; }

   private:
    Variables gvals_;
  };

// -------------------------------------------------------------------------------------------------
}  // namespace ufo
