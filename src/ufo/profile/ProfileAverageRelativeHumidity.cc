/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <set>
#include <sstream>
#include <string>
#include <utility>

#include "ufo/profile/ProfileAverageRelativeHumidity.h"
#include "ufo/profile/ProfileDataHandler.h"
#include "ufo/profile/ProfileDataHolder.h"
#include "ufo/profile/ProfileVerticalAveraging.h"

#include "ufo/utils/metoffice/MetOfficeObservationIDs.h"

namespace ufo {

  static ProfileCheckMaker<ProfileAverageRelativeHumidity>
  makerProfileAverageRelativeHumidity_("AverageRelativeHumidity");

  ProfileAverageRelativeHumidity::ProfileAverageRelativeHumidity
  (const ProfileConsistencyCheckParameters &options)
    : ProfileCheckBase(options)
  {}

  void ProfileAverageRelativeHumidity::runCheck(ProfileDataHandler &profileDataHandler)
  {
    oops::Log::debug() << " Relative humidity averaging" << std::endl;

    // Produce vector of profiles containing data for the wind speed averaging.
    std::vector <std::string> variableNamesInt =
      {ufo::VariableNames::qcflags_relative_humidity,
       ufo::VariableNames::counter_NumGapsRH,
       ufo::VariableNames::extended_obs_space,
       ufo::VariableNames::modellevels_average_relative_humidity_qcflags,
       ufo::VariableNames::InstrType};
    std::vector <std::string> variableNamesFloat =
      {ufo::VariableNames::obs_relative_humidity,
       ufo::VariableNames::pgebd_relative_humidity,
       ufo::VariableNames::LogP_derived,
       ufo::VariableNames::bigPgaps_derived,
       ufo::VariableNames::modellevels_logP_derived,
       ufo::VariableNames::modellevels_logP_rho_derived,
       ufo::VariableNames::modellevels_average_air_temperature_derived,
       ufo::VariableNames::modellevels_average_relative_humidity_derived};
    std::vector <std::string> variableNamesGeoVaLs =
      {ufo::VariableNames::geovals_relative_humidity};

    if (options_.compareWithOPS.value()) {
      variableNamesInt.insert
        (variableNamesInt.end(),
         {"OPS_" +
             std::string(ufo::VariableNames::modellevels_average_relative_humidity_qcflags)});
      variableNamesFloat.insert
        (variableNamesFloat.end(),
         {"OPS_" +
             std::string(ufo::VariableNames::modellevels_average_relative_humidity_derived)});
      variableNamesGeoVaLs.insert
        (variableNamesGeoVaLs.end(),
         {ufo::VariableNames::geovals_average_relative_humidity,
             ufo::VariableNames::geovals_average_relative_humidity_qcflags});
    }

    std::vector <ProfileDataHolder> profiles =
      profileDataHandler.produceProfileVector
      (variableNamesInt,
       variableNamesFloat,
       {},
       variableNamesGeoVaLs,
       {});

    // Run relative humidity averaging on each profile in the original ObsSpace,
    // saving averaged output to the equivalent extended profile.
    const size_t halfnprofs = profileDataHandler.getObsdb().nrecs() / 2;
    for (size_t jprof = 0; jprof < halfnprofs; ++jprof) {
      oops::Log::debug() << "  Profile " << (jprof + 1) << " / " << halfnprofs << std::endl;
      auto& profileOriginal = profiles[jprof];
      auto& profileExtended = profiles[jprof + halfnprofs];
      runCheckOnProfiles(profileOriginal, profileExtended);
    }

    // Fill validation information if required.
    if (options_.compareWithOPS.value()) {
      oops::Log::debug() << " Filling validation data" << std::endl;
      for (size_t jprof = 0; jprof < halfnprofs * 2; ++jprof) {
        fillValidationData(profiles[jprof]);
      }
    }

    // Update data handler with profile information.
    oops::Log::debug() << " Updating data handler" << std::endl;
    profileDataHandler.updateAllProfiles(profiles);
  }

  void ProfileAverageRelativeHumidity::runCheckOnProfiles(ProfileDataHolder &profileOriginal,
                                                          ProfileDataHolder &profileExtended)
  {
    // Check the two profiles are in the correct section of the ObsSpace.
    profileOriginal.checkObsSpaceSection(ufo::ObsSpaceSection::Original);
    profileExtended.checkObsSpaceSection(ufo::ObsSpaceSection::Extended);

    const size_t numProfileLevels = profileOriginal.getNumProfileLevels();
    // Do not perform averaging if there is just one reported level.
    if (numProfileLevels <= 1)
      return;

    const std::vector <float> &rhObs =
      profileOriginal.get<float>(ufo::VariableNames::obs_relative_humidity);
    std::vector <float> &rhPGEBd =
      profileOriginal.get<float>(ufo::VariableNames::pgebd_relative_humidity);
    std::vector <int> &rhFlags =
      profileOriginal.get<int>(ufo::VariableNames::qcflags_relative_humidity);
    // Optional: vector of sonde instrument types.
    const std::vector <int> &InstrType =
      profileOriginal.get<int>(ufo::VariableNames::InstrType);
    std::vector <int> &NumGapsRH =
      profileOriginal.get<int>(ufo::VariableNames::counter_NumGapsRH);

    if (!oops::allVectorsSameNonZeroSize(rhObs, rhPGEBd, rhFlags)) {
      oops::Log::warning() << "At least one vector is the wrong size. "
                           << "Relative humidity averaging will not be performed." << std::endl;
      oops::Log::warning() << "Vector sizes: "
                           << oops::listOfVectorSizes(rhObs,
                                                      rhPGEBd,
                                                      rhFlags)
                           << std::endl;
      // todo(ctgh): Revisit this (and other routines in which a similar choice has been made)
      // when the organisation of the input data becomes clearer.
      return;
    }

    // Obtain GeoVaLs.
    std::vector <float> &geovals_relative_humidity =
      profileOriginal.getGeoVaLVector(ufo::VariableNames::geovals_relative_humidity);
    if (geovals_relative_humidity.empty())
      throw eckit::BadValue("GeoVaLs vector is empty.", Here());

    // Obtain vectors that were produced in the AveragePressure routine.
    const std::vector <float> &LogPA =
      profileExtended.get<float>(ufo::VariableNames::modellevels_logP_rho_derived);
    const std::vector <float> &LogPB =
      profileExtended.get<float>(ufo::VariableNames::modellevels_logP_derived);
    const std::vector <float> &RepLogP =
      profileOriginal.get<float>(ufo::VariableNames::LogP_derived);
    const std::vector <float> &BigGap =
      profileOriginal.get<float>(ufo::VariableNames::bigPgaps_derived);
    if (LogPA.empty() || LogPB.empty() ||
        !oops::allVectorsSameNonZeroSize(RepLogP, BigGap)) {
      std::stringstream errorMessage;
      errorMessage << "At least one model-level vector is the wrong size. "
                   << "Ensure that the AveragePressure routine has been run before this one."
                   << std::endl;
      errorMessage << "Vector sizes: "
                   << oops::listOfVectorSizes(LogPA, LogPB,
                                              RepLogP, BigGap)
                   << std::endl;
      throw eckit::BadValue(errorMessage.str(), Here());
    }

    // Flag reported value if the probability of gross error is too large.
    // Values which have been flagged here, or previously, are not used in the averaging routines.
    for (size_t jlev = 0; jlev < numProfileLevels; ++jlev)
      if (rhPGEBd[jlev] > options_.AvgRH_PGEskip.value())
        rhFlags[jlev] |= ufo::MetOfficeQCFlags::Elem::FinalRejectFlag;

    // Interpolate Sonde RH onto model levels? The alternative is to perform a vertical average.
    const bool RHinterp = options_.AvgRH_Interp;

    // Select model values of log(pressure) to use depending on whether averaging or
    // interpolation is required.
    const std::vector <float> &LogPmodel = RHinterp ? LogPB : LogPA;

    // Average observed relative humidities onto model levels.
    int NumGaps = 0;  // Number of large gaps in reported profile.
    std::vector <float> rhModObs;
    std::vector <int> rhFlagsModObs;
    // Minimum fraction of a model layer that must have been covered (in the vertical coordinate)
    // by observed values in order for averaging onto that layer to be performed.
    const float SondeDZFraction = options_.AvgRH_SondeDZFraction.value();
    calculateVerticalAverage(rhFlags,
                             rhObs,
                             RepLogP,
                             BigGap,
                             LogPmodel,
                             SondeDZFraction,
                             RHinterp ?
                             ProfileAveraging::Method::Interpolation :
                             ProfileAveraging::Method::Averaging,
                             rhFlagsModObs,
                             rhModObs,
                             NumGaps);

    // Increment relative humidity gap counter if necessary.
    if (NumGaps > 0) NumGapsRH[0]++;

    // If the model values of RH are missing, set the equivalent observed values to missing.
    for (size_t mlev = 0; mlev < geovals_relative_humidity.size(); ++mlev)
      if (geovals_relative_humidity[mlev] == missingValueFloat)
        rhModObs[mlev] = missingValueFloat;

    // If the temperature averaging was previously performed,
    // potentially reject further relative humidity observations.
    const std::vector <float> &tModObs =
      profileExtended.get<float>
      (ufo::VariableNames::modellevels_average_air_temperature_derived);
    if (!tModObs.empty()) {
      // Reject the average relative humidity if the equivalent averaged temperature
      // is less than a particular threshold.
      // There is a default value for the temperature threshold which can be overridden
      // by instrument-dependent values.
      float SondeRHminT = options_.AvgRH_AvgTThreshold.value();
      // Overwrite the value if the instrument is a particular type.
      if (!InstrType.empty()) {
        const int profileInstr = InstrType[0];
        const auto &InstrTThresholds = options_.AvgRH_InstrTThresholds.value();
        if (InstrTThresholds.find(profileInstr) != InstrTThresholds.end())
          SondeRHminT = InstrTThresholds.at(profileInstr);
      }

      // Reject sonde RH?
      bool RejectRH = false;
      for (int mlev = 0; mlev < tModObs.size(); ++mlev) {
        // Once RejectRH is true for one level it is true for the
        // remaining levels in the profile. This is the OPS behaviour.
        if (tModObs[mlev] != missingValueFloat &&
            tModObs[mlev] <= ufo::Constants::t0c + SondeRHminT)
          RejectRH = true;
        if (RejectRH)
          rhFlagsModObs[mlev] |= ufo::MetOfficeQCFlags::Elem::PermRejectFlag;
      }
    }

    // Store the relative humidity averaged onto model levels.
    profileExtended.set<float>
      (ufo::VariableNames::modellevels_average_relative_humidity_derived, std::move(rhModObs));

    // Store the QC flags associated with the relative humidity averaging.
    profileExtended.set<int>
      (ufo::VariableNames::modellevels_average_relative_humidity_qcflags, std::move(rhFlagsModObs));
  }

  void ProfileAverageRelativeHumidity::fillValidationData(ProfileDataHolder &profile)
  {
    // Retrieve, then save, the OPS versions of relative humidity averaged onto model levels,
    // and the QC flags associated with the averaging process.
    profile.set<float>
      ("OPS_" + std::string(ufo::VariableNames::modellevels_average_relative_humidity_derived),
       std::move(profile.getGeoVaLVector
                 (ufo::VariableNames::geovals_average_relative_humidity)));

    // The QC flags are stored as floats but are converted to integers here.
    // Due to the loss of precision, 5 must be added to the missing value.
    const std::vector <float>& average_relative_humidity_qcflags_float =
      profile.getGeoVaLVector
      (ufo::VariableNames::geovals_average_relative_humidity_qcflags);
    std::vector <int> average_relative_humidity_qcflags_int
      (average_relative_humidity_qcflags_float.begin(),
       average_relative_humidity_qcflags_float.end());
    std::replace(average_relative_humidity_qcflags_int.begin(),
                 average_relative_humidity_qcflags_int.end(),
                 -2147483648L,
                 -2147483643L);
    profile.set<int>
      ("OPS_" + std::string(ufo::VariableNames::modellevels_average_relative_humidity_qcflags),
       std::move(average_relative_humidity_qcflags_int));
  }
}  // namespace ufo
