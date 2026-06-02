/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/variabletransforms/Cal_RemapScanPosition.h"

namespace ufo {

/************************************************************************************/
//  Cal_RemapScanPosition
/************************************************************************************/

static TransformMaker<Cal_RemapScanPosition>
    makerCal_RemapScanPosition_("RemapScanPosition");

Cal_RemapScanPosition::Cal_RemapScanPosition(
    const Parameters_ &options,
    const ObsFilterData &data,
    const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
    const std::shared_ptr<ioda::ObsDataVector<float>> &obserr)
    : TransformBase(options, data, flags, obserr), parameters_(options) {}

/************************************************************************************/

void Cal_RemapScanPosition::runTransform(const std::vector<bool> &apply) {
  oops::Log::trace() << " --> Renumber satellite scan position"
            << std::endl;
  oops::Log::trace() << "      --> method: " << options_.Method.value() << std::endl;
  oops::Log::trace() << "      --> obsName: " << obsName() << std::endl;

  const size_t nlocs = obsdb_.nlocs();

  std::vector<int> original_scan_position;
  getObservation("MetaData", "sensorScanPosition", original_scan_position, true);

  std::vector<int> remapped_scan_position(nlocs);
  remapped_scan_position.assign(nlocs, missingValueInt);

  // Loop over all obs
  for (size_t jobs = 0; jobs < nlocs; ++jobs) {
    // if the data have been excluded by the where statement
    if (!apply[jobs]) continue;

    // Calculate revised scan position
    if (original_scan_position[jobs] != missingValueInt) {
       remapped_scan_position[jobs] = formulas::RenumberScanPosition(
         original_scan_position[jobs], parameters_.numFOV.value(), parameters_.floorRemap.value());
    }
  }
  // Overwrite variable at existing locations
  obsdb_.put_db("MetaData", "sensorScanPosition", remapped_scan_position);
}
}  // namespace ufo

