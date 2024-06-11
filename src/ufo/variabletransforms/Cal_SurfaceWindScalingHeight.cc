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

#include "ufo/GeoVaLs.h"
#include "ufo/variabletransforms/Cal_SurfaceWindScalingHeight.h"

namespace ufo {

// -------------------------------------------------------------------------------------------------

  // Add class to the factory
  static TransformMaker<Cal_SurfaceWindScalingHeight>
         makerCal_SurfaceWindScalingHeight_("SurfaceWindScalingHeight");

// -------------------------------------------------------------------------------------------------

  Cal_SurfaceWindScalingHeight::Cal_SurfaceWindScalingHeight(
                                const Parameters_ &options,
                                const ObsFilterData &data,
                                const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
                                const std::shared_ptr<ioda::ObsDataVector<float>> &obserr) :
    TransformBase(options, data, flags, obserr), gvals_(),
    heightVariableGroup_(options.heightVariableGroup),
    heightVariableName_(options.heightVariableName)
  {
    oops::Log::trace() << "Cal_SurfaceWindScalingHeight::Constructor start" << std::endl;
    // List of GeoVaLs this transform will need access to
    gvals_ += Variable("GeoVaLs/geopotential_height");
    gvals_ += Variable("GeoVaLs/wind_reduction_factor_at_10m");
    oops::Log::trace() << "Cal_SurfaceWindScalingHeight::Constructor done" << std::endl;
  }

// -------------------------------------------------------------------------------------------------

  void Cal_SurfaceWindScalingHeight::runTransform(const std::vector<bool> &apply) {
    oops::Log::trace() << "Cal_SurfaceWindScalingHeight::runTransform start" << std::endl;

    // Number of observation locations
    const size_t nlocs = obsdb_.nlocs();

    // Initialize the output variable
    std::vector<double> surfaceWindScalingHeight(nlocs);

    // If no locations to process, add to obs space and return
    if (nlocs == 0) {
      obsdb_.put_db("DerivedVariables", "SurfaceWindScalingHeight", surfaceWindScalingHeight);
      return;
    }

    // Pointer to the GeoVaLs
    const ufo::GeoVaLs * gvals = data_.getGeoVaLs();

    // Number of 'full' model levels
    size_t nlevs = gvals->nlevs(oops::Variable{"geopotential_height"});


    // Get observation height and pressure
    std::vector<double> heightVariable(nlocs), stationElevation(nlocs);
    getObservation(heightVariableGroup_, heightVariableName_, heightVariable, true);
    getObservation("MetaData", "stationElevation", stationElevation, true);

    // Missing values
    const double missing = util::missingValue<double>();

    // Containers for GeoVaLs profiles
    std::vector<double> windReductionFactorAt10m(1), surfacePressure(1), surfaceGeometricHeight(1),
                        geopotentialHeight(nlevs);

    // Loop over locations and compute the surface wind scaling factor for height coordinate
    // -------------------------------------------------------------------------------------
    for (int iloc = 0; iloc < nlocs; ++iloc) {
      // Get the GeoVaLs at this location
      gvals->getAtLocation(windReductionFactorAt10m, oops::Variable{"wind_reduction_factor_at_10m"},
                                                                    iloc);
      gvals->getAtLocation(geopotentialHeight, oops::Variable{"geopotential_height"}, iloc);

      // For values above lowest model level the scaling factor is 1
      surfaceWindScalingHeight[iloc] = 1.0;

      // Set to missing if observation is not present
      if (heightVariable[iloc] == missing) {
        continue;
      }

      // Compute the wind scaling factor for height coordinate
      if (heightVariable[iloc] <= geopotentialHeight[nlevs-1]) {
        surfaceWindScalingHeight[iloc] = windReductionFactorAt10m[0];

        if (heightVariable[iloc] < 10.0) {
          surfaceWindScalingHeight[iloc] = windReductionFactorAt10m[0] *
                                           std::max(heightVariable[iloc], 0.0) / 10.0;
        } else if (heightVariable[iloc] > 10.0) {
          surfaceWindScalingHeight[iloc] = 1.0 + (windReductionFactorAt10m[0] - 1.0) *
                                           (geopotentialHeight[nlevs-1] - heightVariable[iloc]) /
                                           (geopotentialHeight[nlevs-1] - 10.0);
        }
      }
    }

    // Add the output to the ObsSpace
    // ------------------------------
    obsdb_.put_db("DerivedVariables", "SurfaceWindScalingHeight", surfaceWindScalingHeight);

    oops::Log::trace() << "Cal_SurfaceWindScalingHeight::runTransform done" << std::endl;
  }

// -------------------------------------------------------------------------------------------------
}  // namespace ufo

