// *
// * (C) Crown copyright 2021, Met Office
// *
// * This software is licensed under the terms of the Apache Licence Version 2.0
// * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
//

#include "ufo/filters/obsfunctions/OrbitAngle.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/Constants.h"

namespace ufo {

static ObsFunctionMaker<OrbitAngle> maker_("OrbitAngle");

OrbitAngle::OrbitAngle(const eckit::LocalConfiguration & conf) {
  // List of required ObsSpace variables
  invars_ += Variable("MetaData/ephemerisLatitude1");
  invars_ += Variable("MetaData/ephemerisLongitude1");
  invars_ += Variable("MetaData/ephemerisLatitude2");
  invars_ += Variable("MetaData/ephemerisLongitude2");
  invars_ += Variable("MetaData/latitude");
  invars_ += Variable("MetaData/longitude");
  invars_ += Variable("MetaData/dateTime");
}

void OrbitAngle::compute(const ObsFilterData & in, ioda::ObsDataVector<float> & out) const {
  // constants
  const float missingFloat = util::missingValue<float>();
  const util::DateTime missingDateTime = util::missingValue<util::DateTime>();

  const util::DateTime firstjan2000(2000, 1, 1, 0, 0, 0);

  const float ecliptic_npole_RA_hours = 18.0;
  const float ecliptic_npole_dec = 66.55 * Constants::deg2rad;

  const double daysPerSecond = 1.0 / 86400.0;
  const double hoursPerSecond = 1.0 / 3600.0;

  const size_t nlocs = in.nlocs();

  // Inputs
  std::vector<float> ephem_lat1(nlocs), ephem_lat2(nlocs);
  std::vector<float> ephem_lon1(nlocs), ephem_lon2(nlocs);
  std::vector<float> view_lat(nlocs), view_lon(nlocs);
  std::vector<util::DateTime> datetimes(nlocs);

  in.get(Variable("MetaData/ephemerisLatitude1"), ephem_lat1);
  in.get(Variable("MetaData/ephemerisLatitude2"), ephem_lat2);
  in.get(Variable("MetaData/ephemerisLongitude1"), ephem_lon1);
  in.get(Variable("MetaData/ephemerisLongitude2"), ephem_lon2);
  in.get(Variable("MetaData/latitude"), view_lat);
  in.get(Variable("MetaData/longitude"), view_lon);
  in.get(Variable("MetaData/dateTime"), datetimes);

  // Output
  std::vector<float> &orbit_angle = out[0];
  std::fill(orbit_angle.begin(), orbit_angle.end(), missingFloat);

  // Statistics
  size_t numMissingEphemLats = 0;
  size_t numMissingEphemLons = 0;
  size_t numMissingLats = 0;
  size_t numMissingLons = 0;
  size_t numMissingDatetimes = 0;
  size_t numFailedAngleCalc  = 0;
  size_t numOutOfRangeDatetimes = 0;
  size_t numInvalidDateTimes = 0;


  for (size_t loc = 0; loc < nlocs; ++loc) {
    // check for useable input data

    // missing time
    if (datetimes[loc] == missingDateTime) {
      ++numMissingDatetimes;
      oops::Log::debug() << "OrbitAngle: missing datetime encountered for ob " << loc
                         << ". Output set to missing data\n";
      continue;
    }

    // missing lat
    if (view_lat[loc] == missingFloat) {
      ++numMissingLats;
      oops::Log::debug() << "OrbitAngle: missing  viewing latitude encountered for ob " << loc
                         << ". Output set to missing data\n";
      continue;
    }

     // missing lon
    if (view_lon[loc] == missingFloat) {
      ++numMissingLons;
      oops::Log::debug() << "OrbitAngle: missing  viewing longitude encountered for ob " << loc
                         << ". Output set to missing data\n";
      continue;
    }

    // missing ephemeris data
    if (ephem_lon1[loc] == missingFloat || ephem_lon2[loc] == missingFloat) {
      ++numMissingEphemLons;
      oops::Log::debug() << "OrbitAngle: missing  ephemeris longitude encountered for ob " << loc
                         << ". Output set to missing data\n";
      continue;
    }

    if (ephem_lat1[loc] == missingFloat || ephem_lat2[loc] == missingFloat) {
      ++numMissingEphemLats;
      oops::Log::debug() << "OrbitAngle: missing  ephemeris latitude encountered for ob " << loc
                         << ". Output set to missing data\n";
      continue;
    }

    const util::DateTime &datetime = datetimes[loc];

    // convert date to constituent parts
    int year, month, day, hour, minute, second;
    datetime.toYYYYMMDDhhmmss(year, month, day, hour, minute, second);

    // equation is only valid for years between 1950 and 2200
    if (year < 1950 || year > 2200) {
      ++numOutOfRangeDatetimes;
      oops::Log::debug() << "OrbitAngle: date is out of valid range for ob " << loc
                         << ". Output set to missing data\n";
      continue;
    }

    // now work out day of year
    util::DateTime startofyear(year, 1, 1, 0, 0, 0);

    double doy = floor((datetime - startofyear).toSeconds() * daysPerSecond)+1;

    // check date range is meaningful
    if (doy < 1 || doy > 366) {
      ++numInvalidDateTimes;
      oops::Log::debug() << "OrbitAngle: day of year is invalid for ob " << loc
                         << ". Output set to missing data\n";
      continue;
    }

    // convert emphermis positions to cartesian
    std::vector<float> evector1 = OrbitAngle::ll_to_xyz(ephem_lat1[loc], ephem_lon1[loc]);
    std::vector<float> evector2 = OrbitAngle::ll_to_xyz(ephem_lat2[loc], ephem_lon2[loc]);

    // now compute position and vel vectors - these define orbit plane
    std::vector<float> posvector(3);
    for (size_t vloc = 0; vloc < posvector.size(); ++vloc) {
      posvector[vloc] = 0.5*(evector1[vloc]+evector2[vloc]);
    }

    std::vector<float> velvector(3);
    for (size_t vloc = 0; vloc < velvector.size(); ++vloc) {
      velvector[vloc] = evector2[vloc]-evector1[vloc];
    }

    // n.b. a simpler version using jedi datetime format would be
    // long jday=floor((datetime - firstjan2000).toSeconds() / secondsPerDay)+1
    // but i have chosen the OPS method for exact comparison

    util::DateTime startofday(year, month, day, 0, 0, 0);

    double hour_of_day = (datetime - startofday).toSeconds()*hoursPerSecond;

    // check hour of day is meaningful
    if (hour_of_day < 0 || hour_of_day >= 24) {
      ++numInvalidDateTimes;
      oops::Log::debug() << "OrbitAngle: hour is invalid for ob " << loc
                         << ". Output set to missing data\n";
      continue;
    }

    double days_since_2001 = (year - 2001) * static_cast<double>(365.25);

    double jday = floor(days_since_2001) + doy +
                  static_cast<double>(364.5) +
                  static_cast<double>(hour_of_day / 24.0);

    // derive local siderial time
    // this algorithm is also used in many science documents on web
    // e.g. for lunar correction
    double local_sideral_time = static_cast<double>(100.46) +
                                0.985647352 * jday +
                                static_cast<double>(hour_of_day * 15);

    float hour_angle = local_sideral_time * Constants::deg2rad -
                       (ecliptic_npole_RA_hours/24.0) * 2.0* M_PI;


    // now derive the solarplane using local sidereal time
    std::vector<float> solarplane_vector(3);
    solarplane_vector[0] =  std::cos(hour_angle) * std::cos(ecliptic_npole_dec);
    solarplane_vector[1] = -std::sin(hour_angle) * std::cos(ecliptic_npole_dec);
    solarplane_vector[2] = std::sin(ecliptic_npole_dec);

    // cross products  which define each plane
    std::vector<float> orbitplane_vector =
                        OrbitAngle::crossprod_xyz(posvector, velvector);
    std::vector<float> refplane_vector =
                        OrbitAngle::crossprod_xyz(solarplane_vector, orbitplane_vector);

    // use dot product to generate angle between ref plane and satellite pos
    float dpr = OrbitAngle::dotproduct(refplane_vector, posvector);

    // Find the angle between the position (p) of the satellite and the reference vector (r)
    // Recall p.r = |p||r|*cos(angle) therefore,  angle = acos(p.r/|p||r|)  thus,  p.r/|p||r| <= 1
    float modp = OrbitAngle::vectormod(posvector);
    float modr = OrbitAngle::vectormod(refplane_vector);

    // use dot product to generate angle between the solar and satellite pos
    float dps = OrbitAngle::dotproduct(solarplane_vector, posvector);

    if (std::abs(dpr) <= modp*modr) {
      orbit_angle[loc] = std::acos(dpr / (modp * modr)) * Constants::rad2deg;

      // resolve sign ambiguity
      if  (dps < 0) {
        orbit_angle[loc] = 360.0 - orbit_angle[loc];
      }
    } else {
        ++numFailedAngleCalc;
    }
  }

  // Notify about "bad" observations
  if (numMissingEphemLats != 0)
    oops::Log::trace() << "OrbitAngle: " << numMissingEphemLats <<
                          " obs had missing ephemeris latitude\n";
  if (numMissingEphemLons != 0)
    oops::Log::trace() << "OrbitAngle: " << numMissingEphemLons <<
                          " obs had missing ephemeris longitude\n";
  if (numMissingLats != 0)
    oops::Log::trace() << "OrbitAngle: " << numMissingLats <<
                          " obs had missing viewing latitude\n";
  if (numMissingLons != 0)
    oops::Log::trace() << "OrbitAngle: " << numMissingLons <<
                          " obs had missing viewing longitude\n";
  if (numMissingDatetimes != 0)
    oops::Log::trace() << "OrbitAngle: " << numMissingDatetimes <<
                          " obs had missing date times \n";
  if (numOutOfRangeDatetimes != 0)
    oops::Log::trace() << "OrbitAngle: " << numOutOfRangeDatetimes <<
                          " obs had out-of-range datetime\n";
  if (numInvalidDateTimes != 0)
    oops::Log::trace() << "OrbitAngle: " << numInvalidDateTimes <<
                          " obs had invalid datetime\n";
  if (numFailedAngleCalc != 0)
    oops::Log::trace() << "OrbitAngle: " << numFailedAngleCalc <<
                          " obs failed orbit angle calculation\n";
  oops::Log::trace().flush();
}

const ufo::Variables & OrbitAngle::requiredVariables() const {
  return invars_;
}

// Return the modulus of a vector
float OrbitAngle::vectormod(const std::vector<float> &a) const {
  float sumsq = 0.0;
  for (size_t vloc = 0; vloc < a.size(); ++vloc) {
    sumsq += a[vloc]*a[vloc];
  }
  return sqrt(sumsq);
}

// Return the lat long location in cartesian coordinates
std::vector<float> OrbitAngle::ll_to_xyz(const float &lat,
                                        const float &lon) const {
  std::vector<float> xyz(3);
  float coslat = std::cos(lat * Constants::deg2rad);
  xyz[0] = std::cos(lon * Constants::deg2rad) * coslat;
  xyz[1] = std::sin(lon * Constants::deg2rad) * coslat;
  xyz[2] = std::sin(lat * Constants::deg2rad);
  return xyz;
}

// Return the cross product of xyz vector
std::vector<float> OrbitAngle::crossprod_xyz(const std::vector<float> &a,
                                            const std::vector<float> &b) const {
  std::vector<float> crossprod(3);
  crossprod[0] = a[1] * b[2] - a[2] * b[1];
  crossprod[1] = a[2] * b[0] - a[0] * b[2];
  crossprod[2] = a[0] * b[1] - a[1] * b[0];
  return crossprod;
}

// Return the dot product of a vector
float OrbitAngle::dotproduct(const std::vector<float> &a,
                            const std::vector<float> &b) const {
  float dotp = 0.0;
  for (size_t vloc = 0; vloc < a.size(); ++vloc) {
    dotp += a[vloc]*b[vloc];
  }
  return dotp;
}


}  // namespace ufo
