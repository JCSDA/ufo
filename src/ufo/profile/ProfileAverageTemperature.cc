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

#include "ufo/profile/ProfileAverageTemperature.h"
#include "ufo/profile/ProfileDataHandler.h"
#include "ufo/profile/ProfileDataHolder.h"
#include "ufo/profile/ProfileVerticalAveraging.h"

#include "ufo/utils/metoffice/MetOfficeObservationIDs.h"

namespace ufo {

  static ProfileCheckMaker<ProfileAverageTemperature>
  makerProfileAverageTemperature_("AverageTemperature");

  ProfileAverageTemperature::ProfileAverageTemperature
  (const ConventionalProfileProcessingParameters &options)
    : ProfileCheckBase(options)
  {}

  void ProfileAverageTemperature::runCheck(ProfileDataHandler &profileDataHandler)
  {
    oops::Log::debug() << " Temperature averaging" << std::endl;

    // Produce vector of profiles containing data for the temperature averaging.
    std::vector <std::string> variableNamesInt =
      {ufo::VariableNames::qcflags_air_temperature,
       ufo::VariableNames::qcflags_relative_humidity,
       ufo::VariableNames::counter_NumGapsT,
       ufo::VariableNames::extended_obs_space};
    std::vector <std::string> variableNamesFloat =
      {ufo::VariableNames::obs_air_temperature,
       ufo::VariableNames::obscorrection_air_temperature,
       ufo::VariableNames::pge_air_temperature,
       ufo::VariableNames::LogP_derived,
       ufo::VariableNames::bigPgaps_derived,
       ufo::VariableNames::modellevels_ExnerP_rho_derived,
       ufo::VariableNames::modellevels_ExnerP_derived,
       ufo::VariableNames::modellevels_logP_rho_derived,
       ufo::VariableNames::modellevels_logP_derived,
       ufo::VariableNames::modellevels_air_temperature_derived,
       ufo::VariableNames::air_temperature_derived};
    oops::Variables variableNamesGeoVaLs{
      {oops::Variable{ufo::VariableNames::geovals_potential_temperature}}};

    if (options_.compareWithOPS.value()) {
      variableNamesInt.insert
        (variableNamesInt.end(),
         {addOPSPrefix(ufo::VariableNames::qcflags_air_temperature)});
      variableNamesFloat.insert
        (variableNamesFloat.end(),
         {addOPSPrefix(ufo::VariableNames::modellevels_air_temperature_derived),
             addOPSPrefix(ufo::VariableNames::air_temperature_derived)});
      variableNamesGeoVaLs.push_back(oops::Variable{ufo::VariableNames::geovals_air_temperature});
      variableNamesGeoVaLs.push_back(oops::Variable
                                     {ufo::VariableNames::geovals_testreference_air_temperature});
      variableNamesGeoVaLs.push_back(oops::Variable
                               {ufo::VariableNames::geovals_testreference_air_temperature_qcflags});
    }

    std::vector <ProfileDataHolder> profiles =
      profileDataHandler.produceProfileVector
      (variableNamesInt,
       variableNamesFloat,
       {},
       variableNamesGeoVaLs);

    // Run temperature averaging on each profile in the original ObsSpace,
    // saving averaged output to the equivalent extended profile.
    const size_t halfnprofs = profileDataHandler.getObsdb().nrecs() / 2;
    for (size_t jprof = 0; jprof < halfnprofs; ++jprof) {
      oops::Log::debug() << "  Profile " << (jprof + 1) << " / " << halfnprofs << std::endl;
      auto& profileOriginal = profiles[jprof];
      auto& profileExtended = profiles[jprof + halfnprofs];
      runCheckOnProfiles(profileOriginal, profileExtended);
    }

    // Modify filter flags according to values of averaged data.
    ProfileAverageUtils::passNonMissingAveragedObservations
      (profileDataHandler,
       profiles,
       "airTemperature",
       ufo::VariableNames::air_temperature_derived);

    // Fill validation information if required.
    if (options_.compareWithOPS.value()) {
      oops::Log::debug() << " Filling validation data" << std::endl;
      for (size_t jprof = 0; jprof < halfnprofs * 2; ++jprof) {
        ProfileAverageUtils::fillValidationData
          (profiles[jprof],
           jprof >= halfnprofs,
           ufo::VariableNames::air_temperature_derived,
           ufo::VariableNames::qcflags_air_temperature,
           oops::Variable{ufo::VariableNames::geovals_testreference_air_temperature},
           oops::Variable{ufo::VariableNames::geovals_testreference_air_temperature_qcflags});
        // Also fill validation values of model air temperature.
        if (jprof >= halfnprofs) {
          auto & profile = profiles[jprof];
          std::vector <float> & model_air_temperature =
            profile.getGeoVaLVector(oops::Variable{ufo::VariableNames::geovals_air_temperature});
          // Ensure all vectors are the correct size to be saved to the ObsSpace.
          const size_t numModelLevels = profile.getNumProfileLevels();
          model_air_temperature.resize(numModelLevels, missingValueFloat);
          profile.set<float>
            (addOPSPrefix(ufo::VariableNames::modellevels_air_temperature_derived),
             std::move(model_air_temperature));
        }
      }
    }

    // Update data handler with profile information.
    oops::Log::debug() << " Updating data handler" << std::endl;
    profileDataHandler.updateAllProfiles(profiles);
  }

  void ProfileAverageTemperature::runCheckOnProfiles(ProfileDataHolder &profileOriginal,
                                                     ProfileDataHolder &profileExtended)
  {
    // Check the two profiles are in the correct section of the ObsSpace.
    profileOriginal.checkObsSpaceSection(ufo::ObsSpaceSection::Original);
    profileExtended.checkObsSpaceSection(ufo::ObsSpaceSection::Extended);

    const size_t numProfileLevels = profileOriginal.getNumProfileLevels();
    const size_t numModelLevels = profileExtended.getNumProfileLevels();

    // Do not perform averaging if there are fewer than two reported levels.
    // Instead, fill the averaged profile vectors with missing values.
    if (numProfileLevels <= 1) {
      ProfileAverageUtils::setProfileMissing<float>(profileExtended,
        {ufo::VariableNames::modellevels_air_temperature_derived,
         ufo::VariableNames::air_temperature_derived});
      ProfileAverageUtils::setProfileMissing<int>(profileExtended,
        {ufo::VariableNames::qcflags_air_temperature});

      // Store the observed temperature in the vector of derived values.
      // The derived values are initially missing, so performing this action
      // ensures that any filters subsequently run on the original ObsSpace
      // will work correctly.
      ProfileAverageUtils::copyProfileValues<float>(profileOriginal,
                                                    ufo::VariableNames::obs_air_temperature,
                                                    ufo::VariableNames::air_temperature_derived);
      return;
    }

    const std::vector <float> &tObs =
      profileOriginal.get<float>(ufo::VariableNames::obs_air_temperature);
    const std::vector <float> &tPGE =
      profileOriginal.get<float>(ufo::VariableNames::pge_air_temperature);
    std::vector <int> &tFlags =
      profileOriginal.get<int>(ufo::VariableNames::qcflags_air_temperature);
    std::vector <int> &rhFlags =
      profileOriginal.get<int>(ufo::VariableNames::qcflags_relative_humidity);
    const std::vector <float> &tObsCorrection =
      profileOriginal.get<float>(ufo::VariableNames::obscorrection_air_temperature);
    std::vector <int> &NumGapsT =
       profileOriginal.get<int>(ufo::VariableNames::counter_NumGapsT);

    if (!oops::allVectorsSameNonZeroSize(tObs, tPGE, tFlags, rhFlags, tObsCorrection)) {
      std::stringstream errorMessage;
      errorMessage << "At least one vector is the wrong size. "
                   << "Temperature averaging will not be performed." << std::endl;
      errorMessage << "Vector sizes: "
                   << oops::listOfVectorSizes(tObs, tPGE, tFlags, rhFlags, tObsCorrection)
                   << std::endl;
      throw eckit::BadValue(errorMessage.str(), Here());
    }

    // Obtain GeoVaLs potential temperature.
    const std::vector <float> &potempGeoVaLs =
      profileOriginal.getGeoVaLVector(oops::Variable
                                      {ufo::VariableNames::geovals_potential_temperature});
    if (potempGeoVaLs.empty())
      throw eckit::BadValue("Potential temperature GeoVaLs vector is empty.", Here());
    const size_t numTLevels = potempGeoVaLs.size();

    // Obtain vectors that were produced in the AveragePressure routine.
    const std::vector <float> &ExnerPA =
      profileExtended.get<float>(ufo::VariableNames::modellevels_ExnerP_rho_derived);
    const std::vector <float> &ExnerPB =
      profileExtended.get<float>(ufo::VariableNames::modellevels_ExnerP_derived);
    const std::vector <float> &LogPA =
      profileExtended.get<float>(ufo::VariableNames::modellevels_logP_rho_derived);
    const std::vector <float> &LogPB =
      profileExtended.get<float>(ufo::VariableNames::modellevels_logP_derived);
    const std::vector <float> &RepLogP =
      profileOriginal.get<float>(ufo::VariableNames::LogP_derived);
    const std::vector <float> &BigGap =
      profileOriginal.get<float>(ufo::VariableNames::bigPgaps_derived);

    if (!oops::allVectorsSameNonZeroSize(ExnerPA, LogPA) ||
        !oops::allVectorsSameNonZeroSize(ExnerPB, LogPB) ||
        !oops::allVectorsSameNonZeroSize(RepLogP, BigGap)) {
      std::stringstream errorMessage;
      errorMessage << "At least one model-level vector is the wrong size. "
                   << "Ensure that the AveragePressure routine has been run before this one."
                   << std::endl;
      errorMessage << "Vector sizes: "
                   << oops::listOfVectorSizes(ExnerPA, LogPA,
                                              ExnerPB, LogPB,
                                              RepLogP, BigGap)
                   << std::endl;
      throw eckit::BadValue(errorMessage.str(), Here());
    }

    // Apply any corrections to observed temperature.
    std::vector <float> tObsFinal;
    correctVector(tObs, tObsCorrection, tObsFinal);
    // Flag reported value if the probability of gross error is too large.
    // Values which have been flagged here, or previously, are not used in the averaging routines.
    for (size_t jlev = 0; jlev < numProfileLevels; ++jlev) {
      if (tPGE[jlev] > options_.AvgT_PGEskip.value()) {
        tFlags[jlev] |= ufo::MetOfficeQCFlags::Elem::FinalRejectFlag;
        // NB the relative humidity flags are modified in this routine and also
        // in the routine that performs RH averaging.
        rhFlags[jlev] |= ufo::MetOfficeQCFlags::Elem::FinalRejectFlag;
      }
    }

    // Average observed temperatures onto model levels.
    int NumGaps = 0;  //  Number of large gaps in reported profile.
    std::vector <float> tModObs;  // T observations averaged onto model levels.
    std::vector <int> tFlagsModObs;  // Flags associated with the averaging procedure.
    std::vector <float> LogP_Min;  // Min log(pressure) used in layer average.
    std::vector <float> LogP_Max;  // Max log(pressure) used in layer average.
    // Minimum fraction of a model layer that must have been covered (in the vertical coordinate)
    // by observed values in order for averaging onto that layer to be performed.
    const float SondeDZFraction = options_.AvgT_SondeDZFraction.value();
    calculateVerticalAverage(tFlags,
                             tObsFinal,
                             RepLogP,
                             BigGap,
                             LogPA,
                             SondeDZFraction,
                             ProfileAveraging::Method::Averaging,
                             tFlagsModObs,
                             tModObs,
                             NumGaps,
                             &LogP_Max,
                             &LogP_Min);

    // Increment temperature gap counter if necessary.
    if (NumGaps > 0) NumGapsT[0]++;

    // Convert model potential temperature to temperature.
    // The model temperature is used in the partial layer corrections below
    // and in subsequent checks on model levels.
    std::vector <float> tModBkg(numModelLevels, missingValueFloat);
    std::copy(potempGeoVaLs.begin(), potempGeoVaLs.end(), tModBkg.begin());
    for (size_t jlev = 0; jlev < potempGeoVaLs.size(); ++jlev) {
      const float ExnerPBLev = ExnerPB[jlev];
      const float potempLev = potempGeoVaLs[jlev];
      if (ExnerPBLev == missingValueFloat || potempLev == missingValueFloat)
        continue;
      tModBkg[jlev] = ExnerPBLev * potempLev;
    }

    // Recalculate average temperature by taking the thickness of the model layers into account.
    // This procedure uses the values of LogP_Max and LogP_Min that were computed in the
    // calculateVerticalAverage routine.
    // This procedure computes a potential temperature O-B increment using linear interpolation
    // of temperature between the layer boundaries.
    // This increment is added to the background value to produce the averaged observation value.
    const double logPref = std::log(ufo::Constants::pref);
    for (int JLev = 0; JLev < numTLevels; ++JLev) {
      if (tModObs[JLev] == missingValueFloat ||
          LogP_Max[JLev] == missingValueFloat ||
          LogP_Min[JLev] == missingValueFloat ||
          LogP_Max[JLev] == LogP_Min[JLev])
        continue;
      // Difference between between the maximum and minimum values of log(pressure)
      // that were obtained when performing the temperature averaging for this layer.
      const double DLogP = LogP_Max[JLev] - LogP_Min[JLev];
      if (JLev < numTLevels - 1) {  // The current level is below the highest model level.
        // Check whether this model layer is less than 99.5% full, in which case partial layer
        // processing is used.
        if (DLogP < 0.995 * (LogPA[JLev] - LogPA[JLev + 1])) {
          // Lower model level used in the processing.
          const int MLev1 = std::max(JLev - 1, 0);
          // Upper model level used in the processing.
          const int MLev2 = std::min(JLev + 1, static_cast<int>(numTLevels) - 1);
          // If any of the required quantities are missing,
          // set the averaged temperature to the missing value.
          if (tModBkg[MLev1] == missingValueFloat ||
              tModBkg[MLev2] == missingValueFloat ||
              tModBkg[JLev] == missingValueFloat ||
              LogPB[MLev2] == missingValueFloat ||
              LogPA[JLev] == missingValueFloat ||
              LogPA[JLev + 1] == missingValueFloat) {
            tModObs[JLev] = missingValueFloat;
          } else {
            // DExner is guaranteed to be nonzero thanks to the requirement
            // that LogP_Max is not equal to LogP_Min.
            const double DExner =  // Difference between Exner pressures.
              std::exp((LogP_Max[JLev] - logPref) * ufo::Constants::rd_over_cp) -
              std::exp((LogP_Min[JLev] - logPref) * ufo::Constants::rd_over_cp);
            // Compute potential temperature.
            const double potemp = ufo::Constants::rd_over_cp * tModObs[JLev] * DLogP / DExner;
            // Model temperature at level JLev.
            const double TLev =
              potempGeoVaLs[JLev] * (ExnerPA[JLev] - ExnerPA[JLev + 1]) /
              (ufo::Constants::rd_over_cp * (LogPA[JLev] - LogPA[JLev + 1]));
            // Log(P) at model layer midpoint in terms of Log(P).
            const double ModLogP_mid = 0.5 * (LogPA[JLev] + LogPA[JLev + 1]);
            // Model temperature gradient.
            const double TGrad = (tModBkg[MLev1] - tModBkg[MLev2]) / (LogPB[MLev1] - LogPB[MLev2]);
            // Model temperature at level P_Max, computed using midpoint and gradient.
            const double TMax = TLev + (LogP_Min[JLev] - ModLogP_mid) * TGrad;
            // Model temperature at level P_Min, computed using midpoint and gradient.
            const double TMin = TLev + (LogP_Max[JLev] - ModLogP_mid) * TGrad;
            // Model potential temperature for P_Max to P_Min.
            const double potempBk =
              ufo::Constants::rd_over_cp * (TMax + TMin) * DLogP / (2.0 * DExner);
            // Temperature increment for partial layer.
            const double Tinc = (potemp - potempBk) * ExnerPB[JLev];
            // Update averaged temperature with increment.
            tModObs[JLev] = tModBkg[JLev] + Tinc;
          }
        } else {
          // This model layer has been fully covered.
          // Determine difference between Exner pressures using model values.
          const double DExner = ExnerPA[JLev] - ExnerPA[JLev + 1];
          // Compute potential temperature.
          const double potemp =
            ufo::Constants::rd_over_cp * tModObs[JLev] * (LogPA[JLev] - LogPA[JLev + 1]) / DExner;
          // Convert potential temperature back to temperature.
          tModObs[JLev] = potemp * ExnerPB[JLev];
        }
      } else {
        // Highest level to be processed.
        const double DExner =
          std::exp((LogP_Max[JLev] - logPref) * ufo::Constants::rd_over_cp) -
          std::exp((LogP_Min[JLev] - logPref) * ufo::Constants::rd_over_cp);
        // Compute potential temperature.
        const double potemp = ufo::Constants::rd_over_cp * tModObs[JLev] * DLogP / DExner;
        // Convert potential temperature back to temperature.
        tModObs[JLev] = potemp * ExnerPB[JLev];
      }
    }

    // Ensure all vectors are the correct size to be saved to the ObsSpace.
    tModObs.resize(numModelLevels, missingValueFloat);
    tFlagsModObs.resize(numModelLevels, 0);

    // Store the model temperature.
    profileExtended.set<float>
      (ufo::VariableNames::modellevels_air_temperature_derived, std::move(tModBkg));

    // Store the temperature averaged onto model levels.
    profileExtended.set<float>
      (ufo::VariableNames::air_temperature_derived, std::move(tModObs));

    // Store the QC flags associated with temperature averaging.
    profileExtended.set<int>
      (ufo::VariableNames::qcflags_air_temperature, std::move(tFlagsModObs));

    // Store the observed temperature in the vector of derived values.
    // The derived values are initially missing, so performing this action
    // ensures that any filters subsequently run on the original ObsSpace
    // will work correctly.
    // Create a copy to avoid moving from a const vector.
    std::vector<float> tObsToSave = tObs;
    profileOriginal.set<float>
      (ufo::VariableNames::air_temperature_derived, std::move(tObsToSave));
  }
}  // namespace ufo
