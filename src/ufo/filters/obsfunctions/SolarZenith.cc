/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/SolarZenith.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/QCflags.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/Constants.h"

namespace ufo {

namespace {

/// Return a vector whose ith element is set to true if and only if observations of all simulated
/// variables at the ith location have been rejected.
std::vector<bool> identifyRejectedObservations(const ObsFilterData &data) {
  const size_t nlocs = data.nlocs();
  const oops::ObsVariables &simulatedVars = data.obsspace().obsvariables();
  std::vector<bool> rejected(nlocs, true);

  std::vector<int> qcflags(nlocs);
  for (size_t iv = 0; iv < simulatedVars.size(); ++iv) {
    data.get(Variable("QCflagsData/" + simulatedVars[iv]), qcflags);
    for (size_t iloc = 0; iloc < nlocs; ++iloc)
      if (qcflags[iloc] == QCflags::pass)
        rejected[iloc] = false;
  }

  return rejected;
}

}  // namespace

static ObsFunctionMaker<SolarZenith> maker_("SolarZenith");

SolarZenith::SolarZenith(const eckit::LocalConfiguration & conf) {
  options_.validateAndDeserialize(conf);

  // List of required ObsSpace variables
  invars_ += Variable("MetaData/latitude");
  invars_ += Variable("MetaData/longitude");
  invars_ += Variable("MetaData/dateTime");
}

void SolarZenith::compute(const ObsFilterData & in, ioda::ObsDataVector<float> & out) const {
  const float missingFloat = util::missingValue<float>();
  const util::DateTime missingDateTime = util::missingValue<util::DateTime>();

  const int secondsPerDay = 60 * 60 * 24;
  const double centuriesPerDay = 1.0 / 36525.0;
  const double hoursPerSecond = 1.0 / 3600.0;
  const double degreesLongitudePerHour = 15.0;
  const double hoursPerDegreeLongitude = 1.0 / degreesLongitudePerHour;
  const double one_over_360 = 1.0 / 360.0;

  const util::DateTime startOfLastDayOf19thCentury(1899, 12, 31, 0, 0, 0);

  const size_t nlocs = in.nlocs();

  // Inputs
  std::vector<float> lats(nlocs), lons(nlocs);
  std::vector<util::DateTime> datetimes(nlocs);
  in.get(Variable("MetaData/latitude"), lats);
  in.get(Variable("MetaData/longitude"), lons);
  in.get(Variable("MetaData/dateTime"), datetimes);

  std::vector<bool> rejected;
  const bool skipRejected = options_.skipRejected;
  if (skipRejected)
    rejected = identifyRejectedObservations(in);

  // Output
  std::vector<float> &zenith = out[0];
  std::fill(zenith.begin(), zenith.end(), missingFloat);

  // Statistics
  size_t numRejected = 0;
  size_t numMissingLats = 0;
  size_t numMissingLons = 0;
  size_t numMissingDatetimes = 0;
  size_t numOutOfRangeLats = 0;
  size_t numOutOfRangeDatetimes = 0;

  // Values dependent on day only (reused for all consecutive datetimes from the same day)
  util::DateTime dayStart = missingDateTime;
  util::DateTime dayEnd = missingDateTime;
  // Equation of time (more details below)
  double eqnt;
  // Sine and cosine of declination
  double sinDecl, cosDecl;

  for (size_t loc = 0; loc < nlocs; ++loc) {
    if (skipRejected && rejected[loc]) {
      ++numRejected;
      oops::Log::debug() << "SolarZenith: ob " << loc << " has already been rejected"
                         << ". Output set to missing data\n";
      continue;
    }
    if (lats[loc] == missingFloat) {
      ++numMissingLats;
      oops::Log::debug() << "SolarZenith: missing latitude encountered for ob " << loc
                         << ". Output set to missing data\n";
      continue;
    }
    if (lons[loc] == missingFloat) {
      ++numMissingLons;
      oops::Log::debug() << "SolarZenith: missing longitude encountered for ob " << loc
                         << ". Output set to missing data\n";
      continue;
    }
    if (datetimes[loc] == missingDateTime) {
      ++numMissingDatetimes;
      oops::Log::debug() << "SolarZenith: missing datetime encountered for ob " << loc
                         << ". Output set to missing data\n";
      continue;
    }
    if (lats[loc] < -90 || lats[loc] > 90) {
      ++numOutOfRangeLats;
      oops::Log::debug() << "SolarZenith: latitude " << lats[loc] << " of ob " << loc
                         << " is out of range. Output set to missing data\n";
      continue;
    }

    const util::DateTime &datetime = datetimes[loc];

    if (datetime < dayStart || datetime >= dayEnd) {
      // Calculate quantities dependent only on the date (not time).
      int year, month, day, hour, minute, second;
      datetime.toYYYYMMDDhhmmss(year, month, day, hour, minute, second);
      if (year > 1950 && year <= 2200) {
        dayStart = util::DateTime(year, month, day, 0, 0, 0);
        dayEnd = dayStart + util::Duration(secondsPerDay);

        // Day since 31 Dec 1899 ("0 Jan 1900")
        // (Used instead of 1 Jan 1900 since 2000 was a leap year.)
        const std::size_t centuryDay =
            (dayStart - startOfLastDayOf19thCentury).toSeconds() / secondsPerDay;
        const double rcd = centuryDay * centuriesPerDay;  // Fraction of days elapsed this century
        const double rcd2 = rcd * rcd;
        double ydeg = (rcd * 36000.769 + 279.697) * one_over_360;
        ydeg = std::fmod(ydeg, 1.0) * 360.0;
        const double yrad = ydeg * Constants::deg2rad;

        // Compute equation of time (in seconds) for this day
        // (No reference for this but it gives the correct answers
        // when compared with table in Norton's Star Atlas.)
        // The linter protests about extra spaces used for alignment, so is disabled.
        eqnt = - (( 93.0 + 14.23 * rcd - 0.0144 * rcd2) * std::sin(yrad))         // NOLINT
               - ((432.5 - 3.71  * rcd - 0.2063 * rcd2) * std::cos(yrad))         // NOLINT
               + ((596.9 - 0.81  * rcd - 0.0096 * rcd2) * std::sin(2.0 * yrad))   // NOLINT
               - ((  1.4 + 0.28  * rcd)                 * std::cos(2.0 * yrad))   // NOLINT
               + ((  3.8 + 0.6   * rcd)                 * std::sin(3.0 * yrad))   // NOLINT
               + (( 19.5 - 0.21  * rcd - 0.0103 * rcd2) * std::cos(3.0 * yrad))   // NOLINT
               - (( 12.8 - 0.03  * rcd)                 * std::sin(4.0 * yrad));  // NOLINT

        // Get solar declination for given day (radians)
        const double sinalp = std::sin((ydeg - eqnt / 240.0) * Constants::deg2rad);
        const double taneqn = 0.43382 - 0.00027 * rcd;
        const double decl = std::atan(taneqn * sinalp);
        eqnt *= hoursPerSecond;  // Convert to hours

        // Sine and cosine of declination
        sinDecl = std::sin(decl);
        cosDecl = std::cos(decl);
      } else {
        ++numOutOfRangeDatetimes;
        oops::Log::debug() << "SolarZenith: date/time " << datetime << " of ob " << loc
                           << "is out of range. Output set to missing data\n";
        continue;
      }
    }

    const double lat = lats[loc];
    const double lon = lons[loc];

    const double latInRadians = lat * Constants::deg2rad;
    const double sinLat = std::sin(latInRadians);
    const double cosLat = std::cos(latInRadians);

    const std::int64_t secondsSinceDayStart = (datetime - dayStart).toSeconds();
    const double hoursSinceDayStart = secondsSinceDayStart * hoursPerSecond;
    const double localSolarTimeInHours = lon * hoursPerDegreeLongitude + eqnt + hoursSinceDayStart;
    // Local hour angle (when longitude is 0, this is the Greenwich hour angle given in the
    // Air Almanac)
    const double hourAngleInRadians =
        (localSolarTimeInHours * degreesLongitudePerHour + 180.0) * Constants::deg2rad;

    const double sinEv = sinDecl * sinLat + cosDecl * cosLat * std::cos(hourAngleInRadians);
    zenith[loc] = (M_PI / 2 - std::asin(sinEv)) * Constants::rad2deg;
  }

  // Notify about "bad" observations
  if (numRejected != 0)
    oops::Log::trace() << "SolarZenith: " << numMissingLats << " obs had already been rejected\n";
  if (numMissingLats != 0)
    oops::Log::trace() << "SolarZenith: " << numMissingLats << " obs had missing latitude\n";
  if (numMissingLons != 0)
    oops::Log::trace() << "SolarZenith: " << numMissingLats << " obs had missing longitude\n";
  if (numOutOfRangeLats != 0)
    oops::Log::trace() << "SolarZenith: " << numOutOfRangeLats
                       << " obs had out-of-range latitude\n";
  if (numOutOfRangeDatetimes != 0)
    oops::Log::trace() << "SolarZenith: " << numOutOfRangeDatetimes
                       << " obs had out-of-range datetime\n";
  oops::Log::trace().flush();
}

const ufo::Variables & SolarZenith::requiredVariables() const {
  return invars_;
}

}  // namespace ufo
