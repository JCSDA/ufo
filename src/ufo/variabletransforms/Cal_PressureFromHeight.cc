/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/variabletransforms/Cal_PressureFromHeight.h"
#include "ufo/utils/Constants.h"

namespace ufo {
/************************************************************************************/
//  Cal_PressureFromHeightForProfile
/************************************************************************************/

static TransformMaker<Cal_PressureFromHeightForProfile>
    makerCal_PressureFromHeightForProfile_("PressureFromHeightForProfile");

Cal_PressureFromHeightForProfile::Cal_PressureFromHeightForProfile(
    const Parameters_ &options, const ObsFilterData &data,
    const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
    const std::shared_ptr<ioda::ObsDataVector<float>> &obserr)
    : TransformBase(options, data, flags, obserr),
      heightCoord_(options.HeightCoord),
      pressureCoord_(options.PressureCoord),
      pressureGroup_(options.PressureGroup)
{}

/************************************************************************************/

void Cal_PressureFromHeightForProfile::runTransform(const std::vector<bool> &apply) {
  oops::Log::trace() << " --> Retrieve Pressure From Height (Profile)"
            << std::endl;
  oops::Log::trace() << "      --> method: " << method() << std::endl;
  oops::Log::trace() << "      --> formulation: " << formulation() << std::endl;
  oops::Log::trace() << "      --> obsName: " << obsName() << std::endl;

  // Get the right method
  switch (method()) {
    case formulas::MethodFormulation::NCAR:
    case formulas::MethodFormulation::NOAA:
    case formulas::MethodFormulation::UKMO:
    default: {
      methodUKMO(apply);
      break;
    }
  }
}

/************************************************************************************/

void Cal_PressureFromHeightForProfile::methodUKMO(const std::vector<bool> &apply) {
  std::vector<float> airTemperature;
  std::vector<float> airTemperatureSurface;
  std::vector<float> geopotentialHeight;
  std::vector<float> dewPointTemperature;
  std::vector<float> dewPointTemperatureSurface;
  std::vector<float> relativeHumidity;
  std::vector<float> relativeHumiditySurface;
  std::vector<float> pressureStation;
  std::vector<float> stationElevation;
  std::vector<float> airPressure;

  float Pvap = missingValueFloat;   // Vapour pressure
  float Zcurrent = missingValueFloat;   // Current height value
  float Tcurrent = missingValueFloat;   // Current temperature value
  float Pprev = missingValueFloat;  // Previous pressure value [ps]
  float Zprev = missingValueFloat;  // Previous height value [m]
  float Tprev = missingValueFloat;  // Previous temperature value [k]

  bool hasBeenUpdated = false;

  const size_t nlocs_ = obsdb_.nlocs();
  ioda::ObsSpace::RecIdxIter irec;

  // Here we can only use data that have not been QCed
  // so making sure UseValidDataOnly_ is set to True
  SetUseValidDataOnly(true);

  // 0. Innitialise the ouput array
  // -------------------------------------------------------------------------------
  getObservation(pressureGroup_, pressureCoord_,
                 airPressure);
  if (airPressure.empty()) {
    airPressure = std::vector<float>(nlocs_);
    std::fill(airPressure.begin(), airPressure.end(), missingValueFloat);
  }

  // 1. get the right variables
  // -------------------------------------------------------------------------------
  // Compulsory meta-data
  getObservation("MetaData", "stationElevation",
                 stationElevation, true);
  // Compulsory surface observation
  getObservation("ObsValue", "stationPressure",
                 pressureStation, true);
  getObservation("ObsValue", "airTemperatureAt2M",
                 airTemperatureSurface, true);

  // Compulsory upper air observation
  getObservation("ObsValue", heightCoord_,
                 geopotentialHeight, true);
  getObservation("ObsValue", "airTemperature",
                 airTemperature, true);

  // Here we have a choice between dew point temperature and relative humidity
  // --> By default we chose dew point temperature first!
  getObservation("ObsValue", "dewPointTemperature",
                 dewPointTemperature);
  if (dewPointTemperature.empty()) {
    // if we don't have dewpoint temperature, use relative humidity.
    getObservation("ObsValue", "relativeHumidity",
                   relativeHumidity, true);
    getObservation("ObsValue", "relativeHumidityAt2M",
                   relativeHumiditySurface, true);
  } else {
    getObservation("ObsValue", "dewPointTemperatureAt2M",
                   dewPointTemperatureSurface, true);
  }

  // 3. making sure we have what we need is here
  // -------------------------------------------------------------------------------
  if (dewPointTemperature.empty()) {
    if (!oops::allVectorsSameNonZeroSize(geopotentialHeight, airTemperature, relativeHumidity)) {
      oops::Log::warning() << "Vector sizes: "
                           << oops::listOfVectorSizes(geopotentialHeight, airTemperature,
                                                      relativeHumidity)
                           << std::endl;
      throw eckit::BadValue("At least one vector is the wrong size or empty out of "
                            "Z, T and Rh", Here());
    }
  } else {
    if (!oops::allVectorsSameNonZeroSize(geopotentialHeight, airTemperature, dewPointTemperature)) {
      oops::Log::warning() << "Vector sizes: "
                           << oops::listOfVectorSizes(geopotentialHeight, airTemperature,
                                                      dewPointTemperature)
                           << std::endl;
      throw eckit::BadValue("At least one vector is the wrong size or empty out of "
                            "Z, T and Td", Here());
    }
  }

  // 4. Starting the calculation
  //    Loop over each record
  // -------------------------------------------------------------------------------------
  for (irec = obsdb_.recidx_begin(); irec != obsdb_.recidx_end(); ++irec) {
    const std::vector<std::size_t> &rSort = obsdb_.recidx_vector(irec);
    size_t ilocs = 0;

    // 4.1 Initialise for surface values
    Pprev = pressureStation[rSort[ilocs]];
    Zprev = stationElevation[rSort[ilocs]];
    Tprev = airTemperatureSurface[rSort[ilocs]];

    // Cycle if stationElevation or airTemperatureSurface or pressureStation is
    // not valid
    if (Zprev == missingValueFloat || Tprev == missingValueFloat ||
        Pprev == missingValueFloat)
      continue;

    // Update Tprev
    if (dewPointTemperature.empty()) {
      // Update Tprev if Rh is valid
      if (relativeHumidity[rSort[ilocs]] != missingValueFloat) {
        Pvap = formulas::SatVaporPres_fromTemp(Tprev, formulation());
        Pvap = formulas::SatVaporPres_correction(Pvap, Tprev, -1.0, formulation());
        Tprev = formulas::VirtualTemp_From_Rh_Psat_P_T(
            relativeHumiditySurface[rSort[ilocs]], Pvap, Pprev, Tprev, formulation());
      }

    } else {
      // Update Tprev if dew point positive
      if (dewPointTemperature[rSort[ilocs]] != missingValueFloat) {
        Pvap = formulas::SatVaporPres_fromTemp(dewPointTemperatureSurface[rSort[ilocs]],
                                      formulation());
        Pvap = formulas::SatVaporPres_correction(Pvap,
                                                 dewPointTemperatureSurface[rSort[ilocs]],
                                                 -1.0,
                                                 formulation());
        Tprev = formulas::VirtualTemp_From_Psat_P_T(Pvap, Pprev, Tprev, formulation());
      }
    }

    // 4.2 Loop over the length of the profile
    for (ilocs = 0; ilocs < rSort.size(); ++ilocs) {
      // if the data have been excluded by the where statement
      if (!apply[rSort[ilocs]]) continue;

      // Take current level values
      Zcurrent = geopotentialHeight[rSort[ilocs]];
      Tcurrent = airTemperature[rSort[ilocs]];

      // Cycle if airPressure is valid
      if (airPressure[rSort[ilocs]] != missingValueFloat) continue;

      // Cycle if geopotentialHeight or airTemperatures is not valid
      if (Zcurrent == missingValueFloat || Tcurrent == missingValueFloat) continue;

      // Update Tcurrent
      if (dewPointTemperature.empty()) {
        // Update Tcurrent if Rh is valid
        if (relativeHumidity[rSort[ilocs]] != missingValueFloat) {
          Pvap = formulas::SatVaporPres_fromTemp(Tprev, formulation());
          Pvap = formulas::SatVaporPres_correction(Pvap, Tprev, -1.0, formulation());
          Tcurrent = formulas::VirtualTemp_From_Rh_Psat_P_T(
              relativeHumidity[rSort[ilocs]], Pvap, Pprev, Tcurrent, formulation());
        }
      } else {
        // Update Tcurrent if dew point positive
        if (dewPointTemperature[rSort[ilocs]] != missingValueFloat) {
          Pvap = formulas::SatVaporPres_fromTemp(dewPointTemperature[rSort[ilocs]],
                                                 formulation());
          Pvap = formulas::SatVaporPres_correction(Pvap,
                                                   dewPointTemperature[rSort[ilocs]],
                                                   -1.0,
                                                   formulation());
          Tcurrent = formulas::VirtualTemp_From_Psat_P_T(Pvap, Pprev, Tcurrent, formulation());
        }
      }
      // Hydrostatic equation:
      airPressure[rSort[ilocs]] =
          Pprev * std::exp((Zprev - Zcurrent) * Constants::grav * 2.0 /
                           (Constants::rd * (Tprev + Tcurrent)));

      // update previous level for next level up
      Pprev = airPressure[rSort[ilocs]];
      Zprev = Zcurrent;
      Tprev = Tcurrent;
      hasBeenUpdated = true;
    }
  }

  if (hasBeenUpdated) {
    // if updated the airPressure
    // assign the derived air pressure as DerivedObsValue
    putObservation(pressureCoord_, airPressure,
                   getDerivedGroup(pressureGroup_));
  }
}

/************************************************************************************/
//  Cal_PressureFromHeightForICAO
/************************************************************************************/
static TransformMaker<Cal_PressureFromHeightForICAO>
    makerCal_PressureFromHeightForICAO_("PressureFromHeightForICAO");

Cal_PressureFromHeightForICAO::Cal_PressureFromHeightForICAO(
    const Parameters_ &options, const ObsFilterData &data,
    const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
    const std::shared_ptr<ioda::ObsDataVector<float>> &obserr)
    : TransformBase(options, data, flags, obserr),
      heightCoord_(options.HeightCoord),
      pressureCoord_(options.PressureCoord),
      pressureGroup_(options.PressureGroup)
{}

/************************************************************************************/

void Cal_PressureFromHeightForICAO::runTransform(const std::vector<bool> &apply) {
  oops::Log::trace() << " Retrieve Pressure From Height (ICAO)" << std::endl;
  oops::Log::trace() << "      --> method: " << method() << std::endl;
  oops::Log::trace() << "      --> formulation: " << formulation() << std::endl;
  oops::Log::trace() << "      --> obsName: " << obsName() << std::endl;

  // Get the right method
  switch (method()) {
    case formulas::MethodFormulation::NCAR:
    case formulas::MethodFormulation::NOAA:
    case formulas::MethodFormulation::UKMO:
    default: {
      methodUKMO(apply);
      break;
    }
  }
}
/************************************************************************************/

void Cal_PressureFromHeightForICAO::methodUKMO(const std::vector<bool> &apply) {
  std::vector<float> geopotentialHeight;
  std::vector<float> airPressure;
  std::vector<float> airPressure_ref;
  bool hasBeenUpdated = false;

  const size_t nlocs_ = obsdb_.nlocs();

  ioda::ObsSpace::RecIdxIter irec;

  // 0. Initialise the ouput array
  // -------------------------------------------------------------------------------
  getObservation(pressureGroup_, pressureCoord_,
                 airPressure);
  if (airPressure.empty()) {
    airPressure = std::vector<float>(nlocs_);
    std::fill(airPressure.begin(), airPressure.end(), missingValueFloat);
  }

  // 1. get the right variables
  // -------------------------------------------------------------------------------
  getObservation("ObsValue", heightCoord_,
                 geopotentialHeight);

  // 2. making sure what we need is here
  // -------------------------------------------------------------------------------
  if (oops::anyVectorEmpty(geopotentialHeight)) {
    oops::Log::warning() << "GeopotentialHeight vector is empty. "
                         << "Check will not be performed." << std::endl;
    throw eckit::BadValue("GeopotentialHeight vector is the wrong size or empty ", Here());
  }

  // 3. Loop over each record
  // -------------------------------------------------------------------------------------
  for (irec = obsdb_.recidx_begin(); irec != obsdb_.recidx_end(); ++irec) {
    const std::vector<std::size_t> &rSort = obsdb_.recidx_vector(irec);
    size_t ilocs = 0;

    // 3.1 Loop over each record
    for (ilocs = 0; ilocs < rSort.size(); ++ilocs) {
      // if the data have been excluded by the where statement
      if (!apply[rSort[ilocs]]) continue;

      // Cycle if airPressure is valid
      if (airPressure[rSort[ilocs]] != missingValueFloat) continue;

      airPressure[rSort[ilocs]] = formulas::Height_To_Pressure_ICAO_atmos(
          geopotentialHeight[rSort[ilocs]], formulation());

      hasBeenUpdated = true;
    }
  }

  if (hasBeenUpdated) {
    // if updated the airPressure
    // assign the derived air pressure as DerivedObsValue
    putObservation(pressureCoord_, airPressure,
                   getDerivedGroup(pressureGroup_));
  }
}
}  // namespace ufo

