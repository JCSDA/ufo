/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/variabletransforms/Cal_Wind.h"
#include "ufo/utils/Constants.h"

#include "ufo/filters/VariableTransformParametersBase.h"


namespace ufo {

/************************************************************************************/
//  Cal_WindSpeedAndDirection
/************************************************************************************/

static TransformMaker<Cal_WindSpeedAndDirection>
    makerCal_WindSpeedAndDirection_("WindSpeedAndDirection");

Cal_WindSpeedAndDirection::Cal_WindSpeedAndDirection(
    const Parameters_ &options,
    const ObsFilterData &data,
    const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
    const std::shared_ptr<ioda::ObsDataVector<float>> &obserr)
    : TransformBase(options, data, flags, obserr),
      group_(options.group),
      eastwardwindvariable_(options.EastwardWindVariable),
      northwardwindvariable_(options.NorthwardWindVariable)
     {}

/************************************************************************************/

void Cal_WindSpeedAndDirection::runTransform(const std::vector<bool> &apply) {
  oops::Log::trace() << " --> Retrieve wind speed and direction"
            << std::endl;
  oops::Log::trace() << "      --> method: " << method() << std::endl;
  oops::Log::trace() << "      --> obsName: " << obsName() << std::endl;

  const size_t nlocs = obsdb_.nlocs();

  std::vector<float> u, v;
  getObservation(group_, eastwardwindvariable_,
                 u, true);
  getObservation(group_, northwardwindvariable_,
                 v, true);

  if (!oops::allVectorsSameNonZeroSize(u, v)) {
    oops::Log::warning() << "Vector sizes: "
                         << oops::listOfVectorSizes(u, v)
                         << std::endl;
    throw eckit::BadValue("At least one vector is the wrong size or empty out of "
                          "U, and V", Here());
  }

  std::vector<float> windSpeed(nlocs), windFromDirection(nlocs);
  windSpeed.assign(nlocs, missingValueFloat);
  windFromDirection.assign(nlocs, missingValueFloat);

  // Loop over all obs
  for (size_t jobs = 0; jobs < nlocs; ++jobs) {
    // if the data have been excluded by the where statement
    if (!apply[jobs]) continue;

    // Calculate wind vector
    if (u[jobs] != missingValueFloat && v[jobs] != missingValueFloat) {
       windFromDirection[jobs] = formulas::GetWindDirection(u[jobs], v[jobs]);
       windSpeed[jobs] = formulas::GetWindSpeed(u[jobs], v[jobs]);
    }
  }
  // put new variable at existing locations
  if (eastwardwindvariable_.find("At10M") != std::string::npos) {
      putObservation("windSpeedAt10M", windSpeed, getDerivedGroup(group_));
      putObservation("windDirectionAt10M", windFromDirection, getDerivedGroup(group_));
  } else {
      putObservation("windSpeed", windSpeed, getDerivedGroup(group_));
      putObservation("windDirection", windFromDirection, getDerivedGroup(group_));
  }
}

/************************************************************************************/
//  Cal_WindComponents
/************************************************************************************/
static TransformMaker<Cal_WindComponents>
    makerCal_WindComponents_("WindComponents");

Cal_WindComponents::Cal_WindComponents(
    const Parameters_ &options,
    const ObsFilterData &data,
    const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
    const std::shared_ptr<ioda::ObsDataVector<float>> &obserr)
  : TransformBase(options, data, flags, obserr),
    group_(options.group),
    windspeedvariable_(options.WindSpeedVariable),
    winddirectionvariable_(options.WindDirectionVariable)
    {}

/************************************************************************************/

void Cal_WindComponents::runTransform(const std::vector<bool> &apply) {
  oops::Log::trace() << " --> Retrieve wind component"
            << std::endl;
  oops::Log::trace() << "      --> method: " << method() << std::endl;
  oops::Log::trace() << "      --> obsName: " << obsName() << std::endl;

  // Note that the dimension label "nchans" is being (mis)used here as a second observed
  // dimension for wind observations. This will be changed in the future.
  const size_t nlocs = obsdb_.nlocs();
  // Parse channels from the config
  const size_t nchans = obsdb_.assimvariables().channels().size();

  if (nchans == 0) {
    std::vector<float> windSpeed, windFromDirection;
    getObservation(group_, windspeedvariable_,
                   windSpeed, true);
    getObservation(group_, winddirectionvariable_,
                   windFromDirection, true);

    if (!oops::allVectorsSameNonZeroSize(windSpeed, windFromDirection)) {
      oops::Log::warning() << "Vector sizes: "
                           << oops::listOfVectorSizes(windSpeed, windFromDirection) << std::endl;
      throw eckit::BadValue("At least one vector is the wrong size or empty out of "
                            "wind speed, and direction", Here());
    }

    std::vector<float> u(nlocs, missingValueFloat), v(nlocs, missingValueFloat);

    for (size_t jobs = 0; jobs < nlocs; ++jobs) {
      // if the data have been excluded by the where statement
      if (!apply[jobs]) continue;

      // Calculate wind vector
      if (windFromDirection[jobs] != missingValueFloat &&
          windSpeed[jobs] != missingValueFloat && windSpeed[jobs] >= 0) {
        u[jobs] = formulas::GetWind_U(windSpeed[jobs], windFromDirection[jobs]);
        v[jobs] = formulas::GetWind_V(windSpeed[jobs], windFromDirection[jobs]);
      }
    }
    if (windspeedvariable_.find("At10M") != std::string::npos) {
      putObservation("windEastwardAt10M", u, getDerivedGroup(group_));
      putObservation("windNorthwardAt10M", v, getDerivedGroup(group_));
    } else {
      putObservation("windEastward", u, getDerivedGroup(group_));
      putObservation("windNorthward", v, getDerivedGroup(group_));
    }
  } else if (nchans > 0) {
    std::vector<int> channels(nchans);
    std::iota(std::begin(channels), std::end(channels), 1);

    std::vector<std::vector<float>> windSpeed(nchans, std::vector<float>(nlocs));
    std::vector<std::vector<float>> windFromDirection(nchans, std::vector<float>(nlocs));
    for (size_t ichan = 0; ichan < nchans; ++ichan) {
      data_.get(Variable(group_ + "/" + windspeedvariable_, channels)[ichan],
                windSpeed[ichan]);
      data_.get(Variable(group_ + "/" + winddirectionvariable_, channels)[ichan],
                windFromDirection[ichan]);
    }

    if (!oops::allVectorsSameNonZeroSize(windSpeed, windFromDirection)) {
      oops::Log::warning() << "Vector sizes: "
                           << oops::listOfVectorSizes(windSpeed, windFromDirection) << std::endl;
      throw eckit::BadValue("At least one vector is the wrong size or empty out of "
                            "wind speed, and direction", Here());
    }

    std::vector<std::vector<float>> u(nchans, std::vector<float>(nlocs, missingValueFloat)),
                                    v(nchans, std::vector<float>(nlocs, missingValueFloat));

    // Loop over all obs
    for (size_t jchan = 0; jchan < nchans; ++jchan) {
      for (size_t jobs = 0; jobs < nlocs; ++jobs) {
        // if the data have been excluded by the where statement
        if (!apply[jobs]) continue;

        // Calculate wind vector
        if (windFromDirection[jchan][jobs] != missingValueFloat &&
              windSpeed[jchan][jobs] != missingValueFloat && windSpeed[jchan][jobs] >= 0) {
          u[jchan][jobs] = formulas::GetWind_U(windSpeed[jchan][jobs],
                           windFromDirection[jchan][jobs]);
          v[jchan][jobs] = formulas::GetWind_V(windSpeed[jchan][jobs],
                           windFromDirection[jchan][jobs]);
        }
      }
    }
    // put new variable at existing locations
    for (size_t jchan = 0; jchan < nchans; ++jchan) {
      std::vector<float> u_chan = u[jchan];
      std::vector<float> v_chan = v[jchan];
      if (windspeedvariable_.find("At10M") != std::string::npos) {
        putObservation("windEastwardAt10M", std::to_string(channels[jchan]), u_chan,
                       {"Location", "Channel"}, getDerivedGroup(group_));
        putObservation("windNorthwardAt10M", std::to_string(channels[jchan]), v_chan,
                       {"Location", "Channel"}, getDerivedGroup(group_));
      } else {
        putObservation("windEastward", std::to_string(channels[jchan]), u_chan,
                       {"Location", "Channel"}, getDerivedGroup(group_));
        putObservation("windNorthward", std::to_string(channels[jchan]), v_chan,
                       {"Location", "Channel"}, getDerivedGroup(group_));
      }
    }
  } else {
    throw eckit::Exception("Channel cannot be negative", Here());
  }
}
}  // namespace ufo
