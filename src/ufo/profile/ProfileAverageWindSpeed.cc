/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>

#include "ufo/profile/ProfileAverageWindSpeed.h"
#include "ufo/profile/ProfileDataHandler.h"
#include "ufo/profile/ProfileDataHolder.h"
#include "ufo/profile/ProfileVerticalAveraging.h"

#include "ufo/utils/metoffice/MetOfficeObservationIDs.h"

namespace ufo {

  static ProfileCheckMaker<ProfileAverageWindSpeed>
  makerProfileAverageWindSpeed_("AverageWindSpeed");

  ProfileAverageWindSpeed::ProfileAverageWindSpeed
  (const ConventionalProfileProcessingParameters &options)
    : ProfileCheckBase(options)
  {}

  void ProfileAverageWindSpeed::runCheck(ProfileDataHandler &profileDataHandler)
  {
    oops::Log::debug() << " Wind speed averaging" << std::endl;

    // Produce vector of profiles containing data for the wind speed averaging.
    std::vector <std::string> variableNamesInt =
      {ufo::ProfileVariableNames::counter_NumGapsU,
       ufo::ProfileVariableNames::counter_NumGapsUWP,
       ufo::ProfileVariableNames::ObsType,
       ufo::ProfileVariableNames::extended_obs_space};
    std::vector <std::string> variableNamesFloat =
      {ufo::ProfileVariableNames::obs_eastward_wind,
       ufo::ProfileVariableNames::obs_northward_wind,
       ufo::ProfileVariableNames::pge_eastward_wind,
       ufo::ProfileVariableNames::LogP_derived,
       ufo::ProfileVariableNames::bigPgaps_derived,
       ufo::ProfileVariableNames::modellevels_logP_derived,
       ufo::ProfileVariableNames::eastward_wind_derived,
       ufo::ProfileVariableNames::northward_wind_derived};
    std::vector <std::string> variableNamesBool =
      {ufo::ProfileVariableNames::diagflags_final_reject_u,
       ufo::ProfileVariableNames::diagflags_perm_reject_u,
       ufo::ProfileVariableNames::diagflags_back_reject_u,
       ufo::ProfileVariableNames::diagflags_interpolation_u,
       ufo::ProfileVariableNames::diagflags_partial_layer_u,
       ufo::ProfileVariableNames::diagflags_final_reject_v,
       ufo::ProfileVariableNames::diagflags_perm_reject_v,
       ufo::ProfileVariableNames::diagflags_back_reject_v,
       ufo::ProfileVariableNames::diagflags_interpolation_v,
       ufo::ProfileVariableNames::diagflags_partial_layer_v};
    oops::Variables variableNamesGeoVaLs{
      {oops::Variable{ufo::ProfileVariableNames::geovals_air_pressure_at_surface}}};

    if (options_.compareWithOPS.value()) {
      variableNamesFloat.insert
        (variableNamesFloat.end(),
         {addOPSPrefix(ufo::ProfileVariableNames::eastward_wind_derived),
             addOPSPrefix(ufo::ProfileVariableNames::northward_wind_derived)});
      variableNamesGeoVaLs.push_back
        (oops::Variable{ufo::ProfileVariableNames::geovals_testreference_eastward_wind});
      variableNamesGeoVaLs.push_back
        (oops::Variable{ufo::ProfileVariableNames::geovals_testreference_eastward_wind_qcflags});
      variableNamesGeoVaLs.push_back
        (oops::Variable{ufo::ProfileVariableNames::geovals_testreference_northward_wind});
      variableNamesGeoVaLs.push_back
        (oops::Variable{ufo::ProfileVariableNames::geovals_testreference_northward_wind_qcflags});
    }

    // In order to correctly handle MPI ranks with zero entries,
    // ensure that all of the variables defined above have been added to the ObsSpace
    // on each rank. This prevents a hang when saving the ObsSpace.
    for (const auto & variableInt : variableNamesInt) {
      const auto & vectorInt = profileDataHandler.get<int>(variableInt);
    }
    for (const auto & variableFloat : variableNamesFloat) {
      const auto & vectorFloat = profileDataHandler.get<float>(variableFloat);
    }
    for (const auto & variableBool : variableNamesBool) {
      const auto & vectorBool = profileDataHandler.get<bool>(variableBool);
    }

    std::vector <ProfileDataHolder> profiles =
      profileDataHandler.produceProfileVector
      (variableNamesInt,
       variableNamesFloat,
       {},
       variableNamesBool,
       variableNamesGeoVaLs);

    // Run wind speed averaging on each profile in the original ObsSpace,
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
       "windEastward",
       ufo::ProfileVariableNames::eastward_wind_derived);
    ProfileAverageUtils::passNonMissingAveragedObservations
      (profileDataHandler,
       profiles,
       "windNorthward",
       ufo::ProfileVariableNames::northward_wind_derived);

    // Fill validation information if required.
    if (options_.compareWithOPS.value()) {
      oops::Log::debug() << " Filling validation data" << std::endl;
      for (size_t jprof = 0; jprof < halfnprofs * 2; ++jprof) {
        ProfileAverageUtils::fillValidationData
          (profiles[jprof],
           jprof >= halfnprofs,
           ufo::ProfileVariableNames::eastward_wind_derived,
           ufo::ProfileVariableNames::diagflags_final_reject_u,
           oops::Variable{ufo::ProfileVariableNames::geovals_testreference_eastward_wind},
           oops::Variable{ufo::ProfileVariableNames::geovals_testreference_eastward_wind_qcflags});
        ProfileAverageUtils::fillValidationData
          (profiles[jprof],
           jprof >= halfnprofs,
           ufo::ProfileVariableNames::northward_wind_derived,
           ufo::ProfileVariableNames::diagflags_final_reject_v,
           oops::Variable{ufo::ProfileVariableNames::geovals_testreference_northward_wind},
           oops::Variable{ufo::ProfileVariableNames::geovals_testreference_northward_wind_qcflags});
      }
    }

    // Update data handler with profile information.
    oops::Log::debug() << " Updating data handler" << std::endl;
    profileDataHandler.updateAllProfiles(profiles);
  }

  void ProfileAverageWindSpeed::runCheckOnProfiles(ProfileDataHolder &profileOriginal,
                                                   ProfileDataHolder &profileExtended)
  {
    // Check the two profiles are in the correct section of the ObsSpace.
    profileOriginal.checkObsSpaceSection(ufo::ObsSpaceSection::Original);
    profileExtended.checkObsSpaceSection(ufo::ObsSpaceSection::Extended);

    const size_t numProfileLevels = profileOriginal.getNumProfileLevels();

    const std::vector<std::string> diagFlagNamesU {
      ufo::ProfileVariableNames::diagflags_final_reject_u,
      ufo::ProfileVariableNames::diagflags_perm_reject_u,
      ufo::ProfileVariableNames::diagflags_back_reject_u,
      ufo::ProfileVariableNames::diagflags_interpolation_u};
    const std::vector<std::string> diagFlagNamesV {
      ufo::ProfileVariableNames::diagflags_final_reject_v,
      ufo::ProfileVariableNames::diagflags_perm_reject_v,
      ufo::ProfileVariableNames::diagflags_back_reject_v,
      ufo::ProfileVariableNames::diagflags_interpolation_v};

    // Do not perform averaging if there are fewer than two reported levels.
    // Instead, fill the averaged profile vectors with missing values.
    if (numProfileLevels <= 1) {
      ProfileAverageUtils::setProfileMissing<float>(profileExtended,
        {ufo::ProfileVariableNames::eastward_wind_derived});
      ProfileAverageUtils::setProfileMissing<float>(profileExtended,
        {ufo::ProfileVariableNames::northward_wind_derived});
      for (const std::string & diagFlagNameU : diagFlagNamesU) {
        ProfileAverageUtils::setProfileTrue(profileExtended, {diagFlagNameU});
      }
      ProfileAverageUtils::setProfileTrue(profileExtended,
        {ufo::ProfileVariableNames::diagflags_partial_layer_u});
      for (const std::string & diagFlagNameV : diagFlagNamesV) {
        ProfileAverageUtils::setProfileTrue(profileExtended, {diagFlagNameV});
      }
      ProfileAverageUtils::setProfileTrue(profileExtended,
        {ufo::ProfileVariableNames::diagflags_partial_layer_v});
      // Store the observed eastward and northward winds in the vectors of derived values.
      // The derived values are initially missing, so performing this action
      // ensures that any filters subsequently run on the original ObsSpace
      // will work correctly.
      ProfileAverageUtils::copyProfileValues<float>
        (profileOriginal,
         ufo::ProfileVariableNames::obs_eastward_wind,
         ufo::ProfileVariableNames::eastward_wind_derived);
      ProfileAverageUtils::copyProfileValues<float>
        (profileOriginal,
         ufo::ProfileVariableNames::obs_northward_wind,
         ufo::ProfileVariableNames::northward_wind_derived);

      return;
    }

    const std::vector <float> &uObs =
      profileOriginal.get<float>(ufo::ProfileVariableNames::obs_eastward_wind);
    const std::vector <float> &vObs =
      profileOriginal.get<float>(ufo::ProfileVariableNames::obs_northward_wind);
    const std::vector <float> &uPGE =
      profileOriginal.get<float>(ufo::ProfileVariableNames::pge_eastward_wind);
    std::vector <int> &NumGapsU =
       profileOriginal.get<int>(ufo::ProfileVariableNames::counter_NumGapsU);
    std::map<std::string, std::vector<bool> > diagFlagVectorsU;
    for (const auto & diagFlagNameU : diagFlagNamesU) {
      const std::vector <bool> diagFlagVector =
        profileOriginal.get<bool>(diagFlagNameU);
      diagFlagVectorsU[diagFlagNameU] = diagFlagVector;
    }
    std::map<std::string, std::vector<bool> > diagFlagVectorsV;
    for (const auto & diagFlagNameV : diagFlagNamesV) {
      const std::vector <bool> diagFlagVector =
        profileOriginal.get<bool>(diagFlagNameV);
      diagFlagVectorsV[diagFlagNameV] = diagFlagVector;
    }
    // Number of gaps for wind profilers.
    std::vector <int> &NumGapsUWP =
       profileOriginal.get<int>(ufo::ProfileVariableNames::counter_NumGapsUWP);
    const std::vector <int> &ObsType =
      profileOriginal.get<int>(ufo::ProfileVariableNames::ObsType);

    if (!oops::allVectorsSameNonZeroSize(uObs, vObs,
                                         uPGE,
                                         ObsType)) {
      std::stringstream errorMessage;
      errorMessage << "At least one vector is the wrong size. "
                   << "Wind speed averaging will not be performed." << std::endl;
      errorMessage << "Vector sizes: "
                   << oops::listOfVectorSizes(uObs, vObs,
                                              uPGE,
                                              ObsType)
                   << std::endl;
      throw eckit::BadValue(errorMessage.str(), Here());
    }

    // Obtain GeoVaLs surface pressure and eastward wind speed.
    std::vector <float> &geovals_air_pressure_at_surface =
      profileOriginal.getGeoVaLVector(oops::Variable{
          ufo::ProfileVariableNames::geovals_air_pressure_at_surface});
    if (geovals_air_pressure_at_surface.empty())
      throw eckit::BadValue("Surface pressure GeoVaLs vector is empty.", Here());

    // Obtain vectors that were produced in the AveragePressure routine.
    const std::vector <float> &LogPB =
      profileExtended.get<float>(ufo::ProfileVariableNames::modellevels_logP_derived);
    const std::vector <float> &RepLogP =
      profileOriginal.get<float>(ufo::ProfileVariableNames::LogP_derived);
    const std::vector <float> &BigGap =
      profileOriginal.get<float>(ufo::ProfileVariableNames::bigPgaps_derived);

    if (LogPB.empty() ||
        !oops::allVectorsSameNonZeroSize(RepLogP, BigGap)) {
      std::stringstream errorMessage;
      errorMessage << "At least one model-level vector is the wrong size. "
                   << "Ensure that the AveragePressure routine has been run before this one."
                   << std::endl;
      errorMessage << "Vector sizes: "
                   << oops::listOfVectorSizes(LogPB,
                                              RepLogP, BigGap)
                   << std::endl;
      throw eckit::BadValue(errorMessage.str(), Here());
    }

    // Create concatenated vector of log(pressure) on both surface and upper-air levels
    // for use in the wind speed averaging.
    std::vector <float> LogPWB = LogPB;
    LogPWB.insert(LogPWB.begin(), std::log(geovals_air_pressure_at_surface[0]));

    // Flag reported value if the probability of gross error is too large.
    // Values which have been flagged here, or previously, are not used in the averaging routines.
    for (size_t jlev = 0; jlev < numProfileLevels; ++jlev) {
      if (uPGE[jlev] > options_.AvgU_PGEskip.value()) {
        diagFlagVectorsU[ufo::ProfileVariableNames::diagflags_final_reject_u][jlev] = true;
        diagFlagVectorsV[ufo::ProfileVariableNames::diagflags_final_reject_v][jlev] = true;
      }
    }

    // Average observed wind speeds onto model levels.
    int NumGaps = 0;  //  Number of large gaps in reported profile
    std::vector <float> uModObs;  // u observations averaged onto model levels.
    // Diagnostic flags associated with the averaging procedure.
    std::map<std::string, std::vector<bool> > diagFlagVectorsUModObs;
    for (const auto & diagFlagNameU : diagFlagNamesU) {
      diagFlagVectorsUModObs[diagFlagNameU] = {};
    }
    std::vector <bool> diagFlagsUPartialLayerModObs;
    // Minimum fraction of a model layer that must have been covered (in the vertical coordinate)
    // by observed values in order for averaging onto that layer to be performed.
    const float SondeDZFraction = options_.AvgU_SondeDZFraction.value();

    calculateVerticalAverage(uObs,
                             RepLogP,
                             BigGap,
                             LogPWB,
                             diagFlagVectorsU,
                             SondeDZFraction,
                             ProfileAveraging::Method::Averaging,
                             uModObs,
                             diagFlagVectorsUModObs,
                             diagFlagsUPartialLayerModObs,
                             NumGaps);

    // Increment wind speed gap counter if necessary.
    if (NumGaps > 0) {
      if (ObsType[0] == ufo::MetOfficeObsIDs::AtmosphericProfile::WindProf)
        NumGapsUWP[0]++;
      else
        NumGapsU[0]++;
    }

    std::vector <float> vModObs;  // v observations averaged onto model levels.
    std::map<std::string, std::vector<bool> > diagFlagVectorsVModObs;
    for (const auto & diagFlagNameV : diagFlagNamesV) {
      diagFlagVectorsVModObs[diagFlagNameV] = {};
    }
    std::vector <bool> diagFlagsVPartialLayerModObs;

    calculateVerticalAverage(vObs,
                             RepLogP,
                             BigGap,
                             LogPWB,
                             diagFlagVectorsV,
                             SondeDZFraction,
                             ProfileAveraging::Method::Averaging,
                             vModObs,
                             diagFlagVectorsVModObs,
                             diagFlagsVPartialLayerModObs,
                             NumGaps);

    // Store the eastward wind speed averaged onto model levels.
    profileExtended.set<float>
      (ufo::ProfileVariableNames::eastward_wind_derived, std::move(uModObs));

    // Store the diagnostic flags associated with the eastward wind averaging.
    for (const auto & diagFlagNameU : diagFlagNamesU) {
      profileExtended.set<bool>
        (diagFlagNameU,
         std::move(diagFlagVectorsUModObs[diagFlagNameU]));
    }
    profileExtended.set<bool>
      (ufo::ProfileVariableNames::diagflags_partial_layer_u,
      std::move(diagFlagsUPartialLayerModObs));

    // Store the northward wind speed averaged onto model levels.
    profileExtended.set<float>
      (ufo::ProfileVariableNames::northward_wind_derived, std::move(vModObs));

    // Store the diagnostic flags associated with the northward wind averaging.
    for (const auto & diagFlagNameV : diagFlagNamesV) {
      profileExtended.set<bool>
        (diagFlagNameV,
         std::move(diagFlagVectorsVModObs[diagFlagNameV]));
    }
    profileExtended.set<bool>
      (ufo::ProfileVariableNames::diagflags_partial_layer_v,
       std::move(diagFlagsVPartialLayerModObs));

    // Store the observed eastward and northward winds in the vectors of derived values.
    // The derived values are initially missing, so performing this action
    // ensures that any filters subsequently run on the original ObsSpace
    // will work correctly.
    // In each case create a copy to avoid moving from a const vector.
    std::vector<float> uObsToSave = uObs;
    profileOriginal.set<float>
      (ufo::ProfileVariableNames::eastward_wind_derived, std::move(uObsToSave));
    std::vector<float> vObsToSave = vObs;
    profileOriginal.set<float>
      (ufo::ProfileVariableNames::northward_wind_derived, std::move(vObsToSave));

    // Save final rejection diagnostic flag vectors, which may have been modified.
    profileOriginal.set<bool>
      (ufo::ProfileVariableNames::diagflags_final_reject_u,
       std::move(diagFlagVectorsU[ufo::ProfileVariableNames::diagflags_final_reject_u]));
    profileOriginal.set<bool>
      (ufo::ProfileVariableNames::diagflags_final_reject_v,
       std::move(diagFlagVectorsV[ufo::ProfileVariableNames::diagflags_final_reject_v]));
  }
}  // namespace ufo
