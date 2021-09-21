/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/variabletransforms/Cal_HeightFromPressure.h"
#include "ufo/utils/Constants.h"

namespace ufo {
/************************************************************************************/
//  Cal_HeightFromPressure
/************************************************************************************/
static TransformMaker<Cal_HeightFromPressure>
    makerCal_HeightFromPressure_("HeightFromPressure");

Cal_HeightFromPressure::Cal_HeightFromPressure(
    const VariableTransformsParameters &options,
    const ObsFilterData &data,
    const std::shared_ptr<ioda::ObsDataVector<int>> &flags)
    : TransformBase(options, data, flags) {}

/************************************************************************************/

void Cal_HeightFromPressure::runTransform(const std::vector<bool> &apply) {
  oops::Log::trace() << " Retrieve Height From Pressure" << std::endl;

  std::vector<float> geopotentialHeight;
  std::vector<float> airPressure;
  bool hasBeenUpdated = false;

  const size_t nlocs_ = obsdb_.nlocs();

  ioda::ObsSpace::RecIdxIter irec;

  // 1. Obtain air pressure from the ObsSpace.
  // Two possible pressure variables are searched for:
  // - air_pressure,
  // - air_pressure_levels.
  // If neither is present an exception is thrown.

  // Do the pressure, and derived height, lie on staggered levels?
  bool staggeredLevels = false;

  getObservation("MetaData", "air_pressure", airPressure);
  if (airPressure.empty()) {
    getObservation("MetaData", "air_pressure_levels", airPressure);
    if (airPressure.empty()) {
      oops::Log::warning() << "Air pressure vector is empty. "
                           << "Check will not be performed." << std::endl;
      throw eckit::BadValue("Air pressure vector is empty ", Here());
    }
    staggeredLevels = true;
  }

  // 2. Initialise the output array
  // -------------------------------------------------------------------------------
  if (staggeredLevels)
    getObservation("MetaData", "geopotential_height_levels", geopotentialHeight);
  else
    getObservation("MetaData", "geopotential_height", geopotentialHeight);

  if (geopotentialHeight.empty()) {
    geopotentialHeight = std::vector<float>(nlocs_);
    std::fill(geopotentialHeight.begin(), geopotentialHeight.end(), missingValueFloat);
  }

  // 3. Loop over each record
  // -------------------------------------------------------------------------------------
  for (irec = obsdb_.recidx_begin(); irec != obsdb_.recidx_end(); ++irec) {
    const std::vector<std::size_t> &rSort = obsdb_.recidx_vector(irec);
    size_t ilocs = 0;

    // 3.1 Loop over each record
    for (ilocs = 0; ilocs < rSort.size(); ++ilocs) {
      // Cycle if the data have been excluded by the where statement
      if (!apply[rSort[ilocs]]) continue;

      // Cycle if geopotential height is valid
      if (geopotentialHeight[rSort[ilocs]] != missingValueFloat) continue;

      geopotentialHeight[rSort[ilocs]] = formulas::Pressure_To_Height(
          airPressure[rSort[ilocs]], method());

      hasBeenUpdated = true;
    }
  }

  if (hasBeenUpdated) {
    // If the geopotential height was updated, save it as a DerivedValue.
    if (staggeredLevels)
      obsdb_.put_db(outputTag, "geopotential_height_levels", geopotentialHeight);
    else
      obsdb_.put_db(outputTag, "geopotential_height", geopotentialHeight);
  }
}
}  // namespace ufo

