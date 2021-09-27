/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/variabletransforms/Cal_ProfileHorizontalDrift.h"
#include "ufo/utils/Constants.h"

namespace ufo {

/************************************************************************************/
//  Cal_ProfileHorizontalDrift
/************************************************************************************/

static TransformMaker<Cal_ProfileHorizontalDrift>
makerCal_ProfileHorizontalDrift_("ProfileHorizontalDrift");

Cal_ProfileHorizontalDrift::Cal_ProfileHorizontalDrift
(const GenericVariableTransformParameters &options,
 const ObsFilterData &data,
 const std::shared_ptr<ioda::ObsDataVector<int>> &flags)
  : TransformBase(options, data, flags) {}

/************************************************************************************/

void Cal_ProfileHorizontalDrift::runTransform(const std::vector<bool> &apply) {
  oops::Log::trace() << " --> Compute horizontal drift lat/lon/time" << std::endl;
  oops::Log::trace() << "      --> method: " << method() << std::endl;

  // Ensure observations have been grouped into profiles.
  if (obsdb_.obs_group_vars().empty())
    throw eckit::UserError("Group variables configuration is empty", Here());

  // Ensure observations have been sorted by air pressure in descending order.
  if (obsdb_.obs_sort_var() != "air_pressure")
    throw eckit::UserError("Sort variable must be air_pressure", Here());
  if (obsdb_.obs_sort_order() != "descending")
    throw eckit::UserError("Profiles must be sorted in descending order", Here());

  // Obtain values from ObsSpace.
  std::vector<float> latitude_in, longitude_in, wind_speed, wind_from_direction, height;
  std::vector<util::DateTime> datetime_in;
  getObservation("MetaData", "latitude", latitude_in, true);
  getObservation("MetaData", "longitude", longitude_in, true);
  getObservation("MetaData", "datetime", datetime_in, true);
  getObservation("ObsValue", "geopotential_height", height, true);
  getObservation("ObsValue", "wind_speed", wind_speed, true);
  getObservation("ObsValue", "wind_from_direction", wind_from_direction, true);

  if (!oops::allVectorsSameNonZeroSize(latitude_in, longitude_in, datetime_in,
                                       height, wind_speed, wind_from_direction)) {
    oops::Log::warning() << "Vector sizes: "
                         << oops::listOfVectorSizes(latitude_in, longitude_in, datetime_in,
                                                    height, wind_speed, wind_from_direction)
                         << std::endl;
    throw eckit::BadValue("At least one vector is the wrong size", Here());
  }

  // Number of locations in the ObsSpace.
  const size_t nlocs = obsdb_.nlocs();

  // Output values are initialised to input values.
  std::vector<float> latitude_out = latitude_in;
  std::vector<float> longitude_out = longitude_in;
  std::vector<util::DateTime> datetime_out = datetime_in;

  // Get correspondence between record numbers and indices in the total sample.
  const std::vector<size_t> &recnums = obsdb_.recidx_all_recnums();

  // Number of profiles in the ObsSpace.
  const size_t nprofs = recnums.size();

  // Perform drift calculation for each profile in the sample.
  for (size_t jprof = 0; jprof < nprofs; ++jprof) {
    const std::vector<size_t> &locs = obsdb_.recidx_vector(recnums[jprof]);
    formulas::horizontalDrift(locs, apply,
                              latitude_in, longitude_in, datetime_in,
                              height, wind_speed, wind_from_direction,
                              latitude_out, longitude_out, datetime_out);
  }

  // Save output values.
  obsdb_.put_db("DerivedMetaData", "latitude", latitude_out);
  obsdb_.put_db("DerivedMetaData", "longitude", longitude_out);
  obsdb_.put_db("DerivedMetaData", "datetime", datetime_out);
}
}  // namespace ufo

