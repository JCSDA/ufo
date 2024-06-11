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
#include "ufo/variabletransforms/Cal_SurfaceWindScalingPressure.h"

namespace ufo {

// -------------------------------------------------------------------------------------------------

  // Add class to the factory
  static TransformMaker<Cal_SurfaceWindScalingPressure>
         makerCal_SurfaceWindScalingPressure_("SurfaceWindScalingPressure");

// -------------------------------------------------------------------------------------------------

  Cal_SurfaceWindScalingPressure::Cal_SurfaceWindScalingPressure(
                                  const GenericVariableTransformParameters &options,
                                  const ObsFilterData &data,
                                  const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
                                  const std::shared_ptr<ioda::ObsDataVector<float>> &obserr) :
  TransformBase(options, data, flags, obserr), gvals_() {
    oops::Log::trace() << "Cal_SurfaceWindScalingPressure::Constructor start" << std::endl;
    // List of GeoVaLs this transform will need access to
    gvals_ += Variable("GeoVaLs/air_pressure");
    gvals_ += Variable("GeoVaLs/wind_reduction_factor_at_10m");
    gvals_ += Variable("GeoVaLs/virtual_temperature");
    gvals_ += Variable("GeoVaLs/surface_pressure");
    oops::Log::trace() << "Cal_SurfaceWindScalingPressure::Constructor done" << std::endl;
  }

// -------------------------------------------------------------------------------------------------

  void Cal_SurfaceWindScalingPressure::runTransform(const std::vector<bool> &apply) {
    oops::Log::trace() << "Cal_SurfaceWindScalingPressure::runTransform start" << std::endl;

    // Number of observation locations
    const size_t nlocs = obsdb_.nlocs();

    // Initialize the output variable
    std::vector<double> surfaceWindScalingPressure(nlocs);

    // If no locations to process, add to obs space and return
    if (nlocs == 0) {
      obsdb_.put_db("DerivedVariables", "SurfaceWindScalingPressure", surfaceWindScalingPressure);
      return;
    }

    // Pointer to the GeoVaLs
    const ufo::GeoVaLs * gvals = data_.getGeoVaLs();

    // Number of 'full' model levels
    size_t nlevs = gvals->nlevs(oops::Variable{"virtual_temperature"});

    // Get observation height and pressure
    std::vector<double> obsPressure(nlocs);
    getObservation("MetaData", "pressure", obsPressure, true);

    // Missing values
    const double missing = util::missingValue<double>();

    // Containers for GeoVaLs profiles
    std::vector<double> windReductionFactorAt10m(1), surfacePressure(1), airPressure(nlevs),
                        virtualTemperature(nlevs);

    // Loop over locations and compute the surface wind scaling factor for pressure coordinate
    // ---------------------------------------------------------------------------------------
    double grav = 9.80665;
    double rd = 2.8705e2;
    double grav_over_rd = grav/rd;

    for (int iloc = 0; iloc < nlocs; ++iloc) {
      // Get the GeoVaLs at this location
      gvals->getAtLocation(windReductionFactorAt10m, oops::Variable{"wind_reduction_factor_at_10m"},
                                                                    iloc);
      gvals->getAtLocation(surfacePressure, oops::Variable{"surface_pressure"}, iloc);
      gvals->getAtLocation(airPressure, oops::Variable{"air_pressure"}, iloc);
      gvals->getAtLocation(virtualTemperature, oops::Variable{"virtual_temperature"}, iloc);

      // For values above lowest model level the scaling factor is 1
      surfaceWindScalingPressure[iloc] = 1.0;

      // Set to missing if observation is not present
      if (obsPressure[iloc] == missing) {
        continue;
      }

      double dLogPressure = std::log(obsPressure[iloc]) - std::log(surfacePressure[0]);
      double rLogPressure = std::log(airPressure[nlevs-1]/surfacePressure[0]);

      if (dLogPressure > rLogPressure) {
        surfaceWindScalingPressure[iloc] = windReductionFactorAt10m[0];
        double dx10 = - grav_over_rd * 10.0 / virtualTemperature[nlevs-1];
        if (dLogPressure < dx10) {
          surfaceWindScalingPressure[iloc] = 1.0 + (windReductionFactorAt10m[0] - 1.0) *
                                             (rLogPressure - dLogPressure) / (rLogPressure - dx10);
        }
      }
    }

    // Add the output to the ObsSpace
    // ------------------------------
    obsdb_.put_db("DerivedVariables", "SurfaceWindScalingPressure", surfaceWindScalingPressure);

    oops::Log::trace() << "Cal_SurfaceWindScalingPressure::runTransform done" << std::endl;
  }

// -------------------------------------------------------------------------------------------------
}  // namespace ufo

