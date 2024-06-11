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
      {ufo::VariableNames::qcflags_eastward_wind,
       ufo::VariableNames::qcflags_northward_wind,
       ufo::VariableNames::counter_NumGapsU,
       ufo::VariableNames::counter_NumGapsUWP,
       ufo::VariableNames::ObsType,
       ufo::VariableNames::extended_obs_space};
    std::vector <std::string> variableNamesFloat =
      {ufo::VariableNames::obs_eastward_wind,
       ufo::VariableNames::obs_northward_wind,
       ufo::VariableNames::pge_eastward_wind,
       ufo::VariableNames::LogP_derived,
       ufo::VariableNames::bigPgaps_derived,
       ufo::VariableNames::modellevels_logP_derived,
       ufo::VariableNames::eastward_wind_derived,
       ufo::VariableNames::northward_wind_derived};
    oops::Variables variableNamesGeoVaLs{
      {oops::Variable{ufo::VariableNames::geovals_surface_pressure}}};

    if (options_.compareWithOPS.value()) {
      variableNamesInt.insert
        (variableNamesInt.end(),
         {addOPSPrefix(ufo::VariableNames::qcflags_eastward_wind),
             addOPSPrefix(ufo::VariableNames::qcflags_northward_wind)});
      variableNamesFloat.insert
        (variableNamesFloat.end(),
         {addOPSPrefix(ufo::VariableNames::eastward_wind_derived),
             addOPSPrefix(ufo::VariableNames::northward_wind_derived)});
      variableNamesGeoVaLs.push_back(oops::Variable
                                     {ufo::VariableNames::geovals_testreference_eastward_wind});
      variableNamesGeoVaLs.push_back(oops::Variable
                                 {ufo::VariableNames::geovals_testreference_eastward_wind_qcflags});
      variableNamesGeoVaLs.push_back(oops::Variable
                                     {ufo::VariableNames::geovals_testreference_northward_wind});
      variableNamesGeoVaLs.push_back(oops::Variable
                                {ufo::VariableNames::geovals_testreference_northward_wind_qcflags});
    }

    std::vector <ProfileDataHolder> profiles =
      profileDataHandler.produceProfileVector
      (variableNamesInt,
       variableNamesFloat,
       {},
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
       ufo::VariableNames::eastward_wind_derived);
    ProfileAverageUtils::passNonMissingAveragedObservations
      (profileDataHandler,
       profiles,
       "windNorthward",
       ufo::VariableNames::northward_wind_derived);

    // Fill validation information if required.
    if (options_.compareWithOPS.value()) {
      oops::Log::debug() << " Filling validation data" << std::endl;
      for (size_t jprof = 0; jprof < halfnprofs * 2; ++jprof) {
        ProfileAverageUtils::fillValidationData
          (profiles[jprof],
           jprof >= halfnprofs,
           ufo::VariableNames::eastward_wind_derived,
           ufo::VariableNames::qcflags_eastward_wind,
           oops::Variable{ufo::VariableNames::geovals_testreference_eastward_wind},
           oops::Variable{ufo::VariableNames::geovals_testreference_eastward_wind_qcflags});
        ProfileAverageUtils::fillValidationData
          (profiles[jprof],
           jprof >= halfnprofs,
           ufo::VariableNames::northward_wind_derived,
           ufo::VariableNames::qcflags_northward_wind,
           oops::Variable{ufo::VariableNames::geovals_testreference_northward_wind},
           oops::Variable{ufo::VariableNames::geovals_testreference_northward_wind_qcflags});
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
    const size_t numModelLevels = profileExtended.getNumProfileLevels();

    // Do not perform averaging if there are fewer than two reported levels.
    // Instead, fill the averaged profile vectors with missing values.
    if (numProfileLevels <= 1) {
      ProfileAverageUtils::setProfileMissing<float>(profileExtended,
        {ufo::VariableNames::eastward_wind_derived});
      ProfileAverageUtils::setProfileMissing<float>(profileExtended,
        {ufo::VariableNames::northward_wind_derived});
      ProfileAverageUtils::setProfileMissing<int>(profileExtended,
        {ufo::VariableNames::qcflags_eastward_wind});
      ProfileAverageUtils::setProfileMissing<int>(profileExtended,
        {ufo::VariableNames::qcflags_northward_wind});

      // Store the observed eastward and northward winds in the vectors of derived values.
      // The derived values are initially missing, so performing this action
      // ensures that any filters subsequently run on the original ObsSpace
      // will work correctly.
      ProfileAverageUtils::copyProfileValues<float>(profileOriginal,
                                                    ufo::VariableNames::obs_eastward_wind,
                                                    ufo::VariableNames::eastward_wind_derived);
      ProfileAverageUtils::copyProfileValues<float>(profileOriginal,
                                                    ufo::VariableNames::obs_northward_wind,
                                                    ufo::VariableNames::northward_wind_derived);

      return;
    }

    const std::vector <float> &uObs =
      profileOriginal.get<float>(ufo::VariableNames::obs_eastward_wind);
    const std::vector <float> &vObs =
      profileOriginal.get<float>(ufo::VariableNames::obs_northward_wind);
    const std::vector <float> &uPGE =
      profileOriginal.get<float>(ufo::VariableNames::pge_eastward_wind);
    std::vector <int> &uFlags =
      profileOriginal.get<int>(ufo::VariableNames::qcflags_eastward_wind);
    std::vector <int> &vFlags =
      profileOriginal.get<int>(ufo::VariableNames::qcflags_northward_wind);
    std::vector <int> &NumGapsU =
       profileOriginal.get<int>(ufo::VariableNames::counter_NumGapsU);
    // Number of gaps for wind profilers.
    std::vector <int> &NumGapsUWP =
       profileOriginal.get<int>(ufo::VariableNames::counter_NumGapsUWP);
    const std::vector <int> &ObsType =
      profileOriginal.get<int>(ufo::VariableNames::ObsType);

    if (!oops::allVectorsSameNonZeroSize(uObs, vObs,
                                         uPGE,
                                         uFlags, vFlags,
                                         ObsType)) {
      std::stringstream errorMessage;
      errorMessage << "At least one vector is the wrong size. "
                   << "Wind speed averaging will not be performed." << std::endl;
      errorMessage << "Vector sizes: "
                   << oops::listOfVectorSizes(uObs, vObs,
                                              uPGE,
                                              uFlags, vFlags,
                                              ObsType)
                   << std::endl;
      throw eckit::BadValue(errorMessage.str(), Here());
    }

    // Obtain GeoVaLs surface pressure and eastward wind speed.
    std::vector <float> &geovals_surface_pressure =
      profileOriginal.getGeoVaLVector(oops::Variable{ufo::VariableNames::geovals_surface_pressure});
    if (geovals_surface_pressure.empty())
      throw eckit::BadValue("Surface pressure GeoVaLs vector is empty.", Here());

    // Obtain vectors that were produced in the AveragePressure routine.
    const std::vector <float> &LogPB =
      profileExtended.get<float>(ufo::VariableNames::modellevels_logP_derived);
    const std::vector <float> &RepLogP =
      profileOriginal.get<float>(ufo::VariableNames::LogP_derived);
    const std::vector <float> &BigGap =
      profileOriginal.get<float>(ufo::VariableNames::bigPgaps_derived);

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
    LogPWB.insert(LogPWB.begin(), std::log(geovals_surface_pressure[0]));

    // Flag reported value if the probability of gross error is too large.
    // Values which have been flagged here, or previously, are not used in the averaging routines.
    for (size_t jlev = 0; jlev < numProfileLevels; ++jlev) {
      if (uPGE[jlev] > options_.AvgU_PGEskip.value()) {
        uFlags[jlev] |= ufo::MetOfficeQCFlags::Elem::FinalRejectFlag;
        vFlags[jlev] |= ufo::MetOfficeQCFlags::Elem::FinalRejectFlag;
      }
    }

    // Average observed wind speeds onto model levels.
    int NumGaps = 0;  //  Number of large gaps in reported profile
    std::vector <float> uModObs;  // u observations averaged onto model levels.
    std::vector <int> uFlagsModObs;  // Flags associated with the u averaging procedure.
    // Minimum fraction of a model layer that must have been covered (in the vertical coordinate)
    // by observed values in order for averaging onto that layer to be performed.
    const float SondeDZFraction = options_.AvgU_SondeDZFraction.value();
    calculateVerticalAverage(uFlags,
                             uObs,
                             RepLogP,
                             BigGap,
                             LogPWB,
                             SondeDZFraction,
                             ProfileAveraging::Method::Averaging,
                             uFlagsModObs,
                             uModObs,
                             NumGaps);

    // Increment wind speed gap counter if necessary.
    if (NumGaps > 0) {
      if (ObsType[0] == ufo::MetOfficeObsIDs::AtmosphericProfile::WindProf)
        NumGapsUWP[0]++;
      else
        NumGapsU[0]++;
    }

    std::vector <float> vModObs;  // v observations averaged onto model levels.
    std::vector <int> vFlagsModObs;  // Flags associated with the v averaging procedure.
    calculateVerticalAverage(vFlags,
                             vObs,
                             RepLogP,
                             BigGap,
                             LogPWB,
                             SondeDZFraction,
                             ProfileAveraging::Method::Averaging,
                             vFlagsModObs,
                             vModObs,
                             NumGaps);

    // Store the eastward wind speed averaged onto model levels.
    profileExtended.set<float>
      (ufo::VariableNames::eastward_wind_derived, std::move(uModObs));

    // Store the QC flags associated with the eastward wind averaging.
    profileExtended.set<int>
      (ufo::VariableNames::qcflags_eastward_wind, std::move(uFlagsModObs));

    // Store the northward wind speed averaged onto model levels.
    profileExtended.set<float>
      (ufo::VariableNames::northward_wind_derived, std::move(vModObs));

    // Store the QC flags associated with the northward wind averaging.
    profileExtended.set<int>
      (ufo::VariableNames::qcflags_northward_wind, std::move(vFlagsModObs));

    // Store the observed eastward and northward winds in the vectors of derived values.
    // The derived values are initially missing, so performing this action
    // ensures that any filters subsequently run on the original ObsSpace
    // will work correctly.
    // In each case create a copy to avoid moving from a const vector.
    std::vector<float> uObsToSave = uObs;
    profileOriginal.set<float>
      (ufo::VariableNames::eastward_wind_derived, std::move(uObsToSave));
    std::vector<float> vObsToSave = vObs;
    profileOriginal.set<float>
      (ufo::VariableNames::northward_wind_derived, std::move(vObsToSave));
  }
}  // namespace ufo
