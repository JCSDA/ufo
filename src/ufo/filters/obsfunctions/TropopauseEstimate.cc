/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/TropopauseEstimate.h"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <sstream>

#include "eckit/exception/Exceptions.h"

#include "ioda/ObsDataVector.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variable.h"

namespace ufo {

static ObsFunctionMaker<TropopauseEstimate> makerObsFuncTropopauseEstimate_("TropopauseEstimate");

// -----------------------------------------------------------------------------

TropopauseEstimate::TropopauseEstimate(const eckit::LocalConfiguration & conf)
  : invars_() {
  oops::Log::debug() << "TropopauseEstimate: config = " << conf << std::endl;
  // Initialize options
  options_.deserialize(conf);

  // We must know the datetime of each observation
  invars_ += Variable("MetaData/dateTime");
  // We must know the latitude of each observation
  invars_ += Variable("MetaData/latitude");

  if (options_.convert_p2z.value())
        oops::Log::debug() << "  TropopauseEstimate: will convert pres to height" << std::endl;
}

// -----------------------------------------------------------------------------

TropopauseEstimate::~TropopauseEstimate() {}

// -----------------------------------------------------------------------------

void TropopauseEstimate::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<float> & out) const {
  const size_t nlocs = in.nlocs();

  // Ensure that only one output variable is expected.
  ASSERT(out.nvars() == 1);

  // Retrieve the options of tropo_equator, tropo_pole, convert_p2z.
  const float tropo_equator = options_.tropo_equator.value();
  const float tropo_pole = options_.tropo_pole.value();
  bool convert_p2z = options_.convert_p2z.value();

  // Retrieve the datetime and latitude.
  std::vector<float> latitude;
  in.get(Variable("MetaData/latitude"), latitude);
  std::vector<util::DateTime> datetimes;
  in.get(Variable("MetaData/dateTime"), datetimes);

  // If datetimes is empty, then we should just exit because there is nothing we can do otherwise.
  if (datetimes.empty()) {
    return;
  }

  int year, month, day, hour, minute, second, day_peak;
  float slope, tropo_start, season_factor;
  std::vector<float> answer(nlocs);

  // Taking a small short-cut.  Rather than loop through all datetimes, we
  // get only the first one to compute the day of the year.
  datetimes[0].toYYYYMMDDhhmmss(year, month, day, hour, minute, second);
  const util::DateTime firstDayOfYear(year, 1, 1, 0, 0, 0);
  const int day_of_year = (datetimes[0] - firstDayOfYear).toSeconds() / (24*3600);

  // Iterate over the observation points.
  for (size_t jj = 0; jj < nlocs; ++jj) {
    if (latitude[jj] > 0.0) {
      day_peak = 203;
    } else {
      day_peak = 21;
    }
    slope = std::max(0.0, (std::abs(latitude[jj])-15.0)/(90.0-15.0));
    tropo_start = (tropo_pole-tropo_equator)*slope + tropo_equator;
    season_factor = std::cos((day_of_year - day_peak)*0.5f*0.0174533f);
    answer[jj] = tropo_start + season_factor*5000.0;

    // If needed, convert pressure to height.
    if (convert_p2z) answer[jj] = 44307.692*(1.0 - pow(answer[jj]/101325.0f, 0.19f));

    out[0][jj] = answer[jj];
  }

  // TODO(gthompsn): Need to add units or other attributes in the future.
  //   varname = "TropopausePressure"; units = "Pa"; or height and meters.
  if (options_.save) {
    out.save("DerivedValue");
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & TropopauseEstimate::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
