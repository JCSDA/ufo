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

class Cal_SurfaceWindScalingHeightParameters: public VariableTransformParametersBase {
  OOPS_CONCRETE_PARAMETERS(Cal_SurfaceWindScalingHeightParameters, VariableTransformParametersBase);

 public:
  oops::Parameter<std::string> heightVariableGroup{"height variable group",
                                                   "DerivedVariables",
                                                   this};
  oops::Parameter<std::string> heightVariableName{"height variable name",
                                                  "adjustedHeight",
                                                  this};
};

// -------------------------------------------------------------------------------------------------

  /*!
  * \brief Calculate the surface wind scaling quantities
  *
  * \details  Within the Variable Transforms filter, apply the transform "SurfaceWindScaling"
  *  in order to calculate the the surface wind scaling for height- and pressure-based coordinates.
  * The formulations are taken from GSI and added to UFO to aid in replicating GSI behavior.
  */

  class Cal_SurfaceWindScalingHeight : public TransformBase {
   public:
    typedef Cal_SurfaceWindScalingHeightParameters Parameters_;
    Cal_SurfaceWindScalingHeight(const Parameters_ &options,
                                 const ObsFilterData &data,
                                 const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
                                 const std::shared_ptr<ioda::ObsDataVector<float>> &obserr);
    // Run variable conversion
    void runTransform(const std::vector<bool> &apply) override;
    Variables requiredVariables() const override { return gvals_; }

   private:
    Variables gvals_;
    std::string heightVariableGroup_;
    std::string heightVariableName_;
  };

// -------------------------------------------------------------------------------------------------
}  // namespace ufo
