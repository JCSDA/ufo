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
#include "ufo/variabletransforms/Cal_AdjustedHeightCoordinate.h"

namespace ufo {

// -------------------------------------------------------------------------------------------------

  // Add class to the factory
  static TransformMaker<Cal_AdjustedHeightCoordinate>
         makerCal_AdjustedHeightCoordinate_("AdjustedHeightCoordinate");

// -------------------------------------------------------------------------------------------------

  Cal_AdjustedHeightCoordinate::Cal_AdjustedHeightCoordinate(
                                const GenericVariableTransformParameters & options,
                                const ObsFilterData &data,
                                const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
                                const std::shared_ptr<ioda::ObsDataVector<float>> &obserr) :
  TransformBase(options, data, flags, obserr), gvals_() {
    oops::Log::trace() << "Cal_AdjustedHeightCoordinate::Constructor start" << std::endl;
    // List of GeoVaLs this transform will need access to
    gvals_ += Variable("GeoVaLs/surface_geometric_height");
    oops::Log::trace() << "Cal_AdjustedHeightCoordinate::Constructor done" << std::endl;
  }

// -------------------------------------------------------------------------------------------------

  void Cal_AdjustedHeightCoordinate::runTransform(const std::vector<bool> &apply) {
    oops::Log::trace() << "Cal_AdjustedHeightCoordinate::runTransform start" << std::endl;

    // Number of observation locations
    const size_t nlocs = obsdb_.nlocs();

    // Initialize the output variable
    std::vector<double> obsHeightNew(nlocs);

    // If no locations to process, add to obs space and return
    if (nlocs == 0) {
      obsdb_.put_db("DerivedVariables", "adjustedHeight", obsHeightNew);
      return;
    }

    // Pointer to the GeoVaLs
    const ufo::GeoVaLs * gvals = data_.getGeoVaLs();

    // Get observation height and pressure
    std::vector<double> obsHeight(nlocs), stationElevation(nlocs);
    getObservation("MetaData", "height", obsHeight, true);
    getObservation("MetaData", "stationElevation", stationElevation, true);

    // Missing values
    const double missing = util::missingValue<double>();

    // Containers for GeoVaLs profiles
    std::vector<double> surfaceGeometricHeight(1);

    // Adjust observation height
    // -------------------------
    for (int iloc = 0; iloc < nlocs; ++iloc) {
      // Set to missing if observation is not present
      if (obsHeight[iloc] == missing) {
        obsHeightNew[iloc] = missing;
        continue;
      }

      // Height above the station
      double heightAboveStation = obsHeight[iloc] - stationElevation[iloc];

      // Set a scaling factor for the difference between the surface height and the station
      // elevation. F is a function that approximately scales between 0 (at the station) and 1
      // (at 1000m above the station)
      double factor;
      if (heightAboveStation > 1000.0) {
        factor = 1.0;
      } else if (heightAboveStation > 10.0) {
        factor = heightAboveStation/990.0;  //- 1.0/99.0 //Should there be a constant -1/99 here?
      } else {
        factor = 0.0;
      }

      // Get the surface geometric height
      gvals->getAtLocation(surfaceGeometricHeight, oops::Variable{"surface_geometric_height"},
                           iloc);

      // Subtract off the scaled station elevation
      // When f is zero (near the station) subtract off the station elevation
      // When f is one  (far above the station) subtract off the surface height
      obsHeightNew[iloc] = obsHeight[iloc] - factor        * surfaceGeometricHeight[0] +
                                            (factor - 1.0) * stationElevation[iloc];
    }

    // Add the output to the ObsSpace
    // ------------------------------
    obsdb_.put_db("DerivedVariables", "adjustedHeight", obsHeightNew);

    oops::Log::trace() << "Cal_AdjustedHeightCoordinate::runTransform done" << std::endl;
  }

// -------------------------------------------------------------------------------------------------
}  // namespace ufo

