/*
 * (C) Crown copyright 2023, Met Office
 * (C) Copyright 2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <cmath>
#include <string>

#include "oops/util/missingValues.h"

#include "ufo/variabletransforms/Cal_SurfaceWindScalingCombined.h"

namespace ufo {

// -------------------------------------------------------------------------------------------------

  // Add class to the factory
  static TransformMaker<Cal_SurfaceWindScalingCombined>
         makerCal_SurfaceWindScalingCombined_("SurfaceWindScalingCombined");

// -------------------------------------------------------------------------------------------------

  Cal_SurfaceWindScalingCombined::Cal_SurfaceWindScalingCombined(
                                const GenericVariableTransformParameters &options,
                                const ObsFilterData &data,
                                const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
                                const std::shared_ptr<ioda::ObsDataVector<float>> &obserr) :
  TransformBase(options, data, flags, obserr), gvals_() {
    oops::Log::trace() << "Cal_SurfaceWindScalingCombined::Constructor start" << std::endl;
    oops::Log::trace() << "Cal_SurfaceWindScalingCombined::Constructor done" << std::endl;
  }

// -------------------------------------------------------------------------------------------------

  void Cal_SurfaceWindScalingCombined::runTransform(const std::vector<bool> &apply) {
    oops::Log::trace() << "Cal_SurfaceWindScalingCombined::runTransform start" << std::endl;

    // Number of observation locations
    const size_t nlocs = obsdb_.nlocs();

    // Initialize the output variable
    std::vector<double> SurfaceWindScalingCombined(nlocs);

    // If no locations to process, add to obs space and return
    if (nlocs == 0) {
      obsdb_.put_db("DerivedVariables", "SurfaceWindScalingCombined", SurfaceWindScalingCombined);
      return;
    }

    // Get observation quantities
    std::vector<double> height(nlocs), SurfaceWindScalingHeight(nlocs),
                        SurfaceWindScalingPressure(nlocs);
    getObservation("MetaData", "height", height, true);
    getObservation("DerivedVariables", "SurfaceWindScalingHeight", SurfaceWindScalingHeight, true);
    getObservation("DerivedVariables", "SurfaceWindScalingPressure", SurfaceWindScalingPressure,
                   true);

    // Missing values
    const double missing = util::missingValue<double>();

    // Loop over locations and compute the surface wind scaling factor for height coordinate
    // -------------------------------------------------------------------------------------
    for (int iloc = 0; iloc < nlocs; ++iloc) {
      if (height[iloc] != missing) {
        // If height is not missing set to height version
        SurfaceWindScalingCombined[iloc] = SurfaceWindScalingHeight[iloc];
      } else {
        // Otherwise set to pressure version
        SurfaceWindScalingCombined[iloc] = SurfaceWindScalingPressure[iloc];
      }
    }

    // Add the output to the ObsSpace
    // ------------------------------
    obsdb_.put_db("DerivedVariables", "SurfaceWindScalingCombined", SurfaceWindScalingCombined);

    oops::Log::trace() << "Cal_SurfaceWindScalingCombined::runTransform done" << std::endl;
  }

// -------------------------------------------------------------------------------------------------
}  // namespace ufo

