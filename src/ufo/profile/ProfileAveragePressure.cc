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

#include "ufo/profile/ProfileAveragePressure.h"

#include "ufo/profile/ProfileDataHandler.h"
#include "ufo/profile/ProfileDataHolder.h"

#include "ufo/utils/metoffice/MetOfficeObservationIDs.h"

namespace ufo {

  static ProfileCheckMaker<ProfileAveragePressure>
  makerProfileAveragePressure_("AveragePressure");

  ProfileAveragePressure::ProfileAveragePressure
  (const ConventionalProfileProcessingParameters &options)
    : ProfileCheckBase(options)
  {}

  void ProfileAveragePressure::logPressure(const std::vector <float> &pressures,
                                           std::vector <float> &logP)
  {
    logP = pressures;
    std::transform(logP.begin(), logP.end(), logP.begin(),
                   [this](float PLev){return PLev > 0 ? std::log(PLev) : missingValueFloat;});
  }

  void ProfileAveragePressure::ExnerPressure(const std::vector <float> &pressures,
                                             std::vector <float> &ExnerP)
  {
    ExnerP = pressures;
    std::transform(ExnerP.begin(), ExnerP.end(), ExnerP.begin(),
                   [this](float PLev){return PLev > 0 ?
                       std::pow(PLev / ufo::Constants::pref, ufo::Constants::rd_over_cp) :
                       missingValueFloat;});
  }

  void ProfileAveragePressure::bigPressureGaps(const std::vector <float> &pressures,
                                               const int ObsType,
                                               std::vector <float> &bigPgaps)
  {
    bigPgaps = pressures;
    const float GapSize = std::log(10.0f);
    const float GapFactor = ObsType == ufo::MetOfficeObsIDs::AtmosphericProfile::WindProf ?
      options_.AvgP_WinProGapFactor.value() * GapSize :
      options_.AvgP_SondeGapFactor.value() * GapSize;
    const float GapLogPDiffMin = options_.AvgP_GapLogPDiffMin.value();
    // Subtracting log(100.0) from PLev effectively converts the pressure to hPa.
    std::transform(bigPgaps.begin(), bigPgaps.end(), bigPgaps.begin(),
                   [this, GapFactor, GapLogPDiffMin](float PLev){return PLev > 0.0 ?
                       GapFactor / std::max(std::log(PLev) - std::log(100.0f), GapLogPDiffMin) :
                       missingValueFloat;});
  }

  void ProfileAveragePressure::runCheck(ProfileDataHandler &profileDataHandler)
  {
    oops::Log::debug() << " Pressure transformation" << std::endl;

    // Produce vector of profiles containing data for the pressure transformation.
    std::vector <std::string> variableNamesFloat =
      {ufo::VariableNames::obs_air_pressure,
       ufo::VariableNames::LogP_derived,
       ufo::VariableNames::bigPgaps_derived,
       ufo::VariableNames::modellevels_logP_rho_derived,
       ufo::VariableNames::modellevels_logP_derived,
       ufo::VariableNames::modellevels_ExnerP_rho_derived,
       ufo::VariableNames::modellevels_ExnerP_derived};
    oops::Variables variableNamesGeoVaLs{{oops::Variable{ufo::VariableNames::geovals_pressure},
       oops::Variable{ufo::VariableNames::geovals_pressure_rho_minus_one}}};

    if (options_.compareWithOPS.value()) {
      variableNamesFloat.insert
        (variableNamesFloat.end(),
         {addOPSPrefix(ufo::VariableNames::modellevels_logP_rho_derived),
             addOPSPrefix(ufo::VariableNames::modellevels_logP_derived),
             addOPSPrefix(ufo::VariableNames::modellevels_ExnerP_rho_derived),
             addOPSPrefix(ufo::VariableNames::modellevels_ExnerP_derived)});
      variableNamesGeoVaLs.push_back(oops::Variable
                                     {ufo::VariableNames::geovals_testreference_logP_rho});
      variableNamesGeoVaLs.push_back(oops::Variable
                                     {ufo::VariableNames::geovals_testreference_logP});
      variableNamesGeoVaLs.push_back(oops::Variable
                                     {ufo::VariableNames::geovals_testreference_ExnerP_rho});
      variableNamesGeoVaLs.push_back(oops::Variable
                                     {ufo::VariableNames::geovals_testreference_ExnerP});
    }

    std::vector <ProfileDataHolder> profiles =
      profileDataHandler.produceProfileVector
      ({ufo::VariableNames::qcflags_observation_report,
          ufo::VariableNames::ObsType,
          ufo::VariableNames::extended_obs_space},
        variableNamesFloat,
        {},
        variableNamesGeoVaLs);

    // Run pressure transformation on each profile in the original ObsSpace,
    // saving transformed output to the equivalent extended profile.
    const size_t halfnprofs = profileDataHandler.getObsdb().nrecs() / 2;
    for (size_t jprof = 0; jprof < halfnprofs; ++jprof) {
      oops::Log::debug() << "  Profile " << (jprof + 1) << " / " << halfnprofs << std::endl;
      auto& profileOriginal = profiles[jprof];
      auto& profileExtended = profiles[jprof + halfnprofs];
      runCheckOnProfiles(profileOriginal, profileExtended);
    }

    // Fill validation information if required.
    if (options_.compareWithOPS.value()) {
      for (size_t jprof = 0; jprof < halfnprofs * 2; ++jprof)
        fillValidationData(profiles[jprof], jprof >= halfnprofs);
    }

    // Update data handler with profile information.
    oops::Log::debug() << " Updating data handler" << std::endl;
    profileDataHandler.updateAllProfiles(profiles);
  }

  void ProfileAveragePressure::runCheckOnProfiles(ProfileDataHolder &profileOriginal,
                                                  ProfileDataHolder &profileExtended)
  {
    // Check the two profiles are in the correct section of the ObsSpace.
    profileOriginal.checkObsSpaceSection(ufo::ObsSpaceSection::Original);
    profileExtended.checkObsSpaceSection(ufo::ObsSpaceSection::Extended);

    const size_t numProfileLevels = profileOriginal.getNumProfileLevels();
    const size_t numModelLevels = profileExtended.getNumProfileLevels();

    const std::vector <float> &pressures =
      profileOriginal.get<float>(ufo::VariableNames::obs_air_pressure);
    const std::vector <int> &ObsType =
      profileOriginal.get<int>(ufo::VariableNames::ObsType);
    const std::vector <int> &ReportFlags =
      profileOriginal.get<int>(ufo::VariableNames::qcflags_observation_report);

    if (!oops::allVectorsSameNonZeroSize(pressures, ObsType, ReportFlags)) {
      std::stringstream errorMessage;
      errorMessage << "At least one vector is the wrong size. "
                   << "Pressure transformation will not be performed." << std::endl;
      errorMessage << "Vector sizes: "
                   << oops::listOfVectorSizes(pressures, ObsType, ReportFlags)
                   << std::endl;
      throw eckit::BadValue(errorMessage.str(), Here());
    }

    // Initialise vectors on reported levels.
    std::vector <float> logP(numProfileLevels, missingValueFloat);  // log(pressure)
    std::vector <float> bigPgaps(numProfileLevels, missingValueFloat);  // Big gaps in pressure

    // Initialise vectors on model rho levels.
    std::vector <float> geovals_testreference_logP_rho
      (numModelLevels, missingValueFloat);  // log(pressure)
    std::vector <float> geovals_testreference_ExnerP_rho
      (numModelLevels, missingValueFloat);  // Exner pressure

    // Initialise vectors on model theta levels.
    std::vector <float> geovals_testreference_logP
      (numModelLevels, missingValueFloat);  // log(pressure)
    std::vector <float> geovals_testreference_ExnerP
      (numModelLevels, missingValueFloat);  // Exner pressure

    // Determine transformed pressures, unless:
    // - there are zero or one reported levels in the profile, or
    // - certain QC flags have previously been set.
    if (numProfileLevels > 1 &&
        !(ReportFlags[0] & ufo::MetOfficeQCFlags::WholeObReport::FinalRejectReport ||
          ReportFlags[0] & ufo::MetOfficeQCFlags::WholeObReport::OutOfAreaReport)) {
      // Determine log(P) on reported levels.
      logPressure(pressures, logP);

      // Determine size of big gaps for each reported level.
      bigPressureGaps(pressures, ObsType[0], bigPgaps);

      // Calculate log(P) and Exner pressure for GeoVaLs on rho levels.
      const std::vector <float> &geovals_pressure_rho =
        profileOriginal.getGeoVaLVector(oops::Variable
                                        {ufo::VariableNames::geovals_pressure_rho_minus_one});
      logPressure(geovals_pressure_rho, geovals_testreference_logP_rho);
      ExnerPressure(geovals_pressure_rho, geovals_testreference_ExnerP_rho);

      // Calculate log(P) and Exner pressure for GeoVaLs on theta levels.
      const std::vector <float> &geovals_pressure =
        profileOriginal.getGeoVaLVector(oops::Variable{ufo::VariableNames::geovals_pressure});
      // Require these GeoVaLs to be present (unlike the case on rho levels).
      if (geovals_pressure.empty())
        throw eckit::BadValue("geovals_pressure is empty", Here());
      logPressure(geovals_pressure, geovals_testreference_logP);
      ExnerPressure(geovals_pressure, geovals_testreference_ExnerP);
    }

    // Store the transformed pressures on reported levels.
    profileOriginal.set<float>(ufo::VariableNames::LogP_derived, std::move(logP));
    profileOriginal.set<float>(ufo::VariableNames::bigPgaps_derived, std::move(bigPgaps));

    // Ensure all vectors are the correct size to be saved to the ObsSpace.
    geovals_testreference_logP.resize(numModelLevels, missingValueFloat);
    geovals_testreference_ExnerP.resize(numModelLevels, missingValueFloat);

    // Store the transformed pressures on model levels.
    profileExtended.set<float>(ufo::VariableNames::modellevels_logP_rho_derived,
                               std::move(geovals_testreference_logP_rho));
    profileExtended.set<float>(ufo::VariableNames::modellevels_logP_derived,
                               std::move(geovals_testreference_logP));
    profileExtended.set<float>(ufo::VariableNames::modellevels_ExnerP_rho_derived,
                               std::move(geovals_testreference_ExnerP_rho));
    profileExtended.set<float>(ufo::VariableNames::modellevels_ExnerP_derived,
                               std::move(geovals_testreference_ExnerP));
  }

  void ProfileAveragePressure::fillValidationData(ProfileDataHolder &profile,
                                                  bool extended_obs_space)
  {
    // Retrieve, then save, the OPS versions of the transformed pressures on model levels.
    // The quantities retrieved depend on whether the profile lies in the original
    // or averaged sections of the ObsSpace.
    if (extended_obs_space) {
      profile.set<float>
        (addOPSPrefix(ufo::VariableNames::modellevels_logP_rho_derived),
         std::move(profile.getGeoVaLVector
                   (oops::Variable{ufo::VariableNames::geovals_testreference_logP_rho})));
      profile.set<float>
        (addOPSPrefix(ufo::VariableNames::modellevels_logP_derived),
         std::move(profile.getGeoVaLVector
                   (oops::Variable{ufo::VariableNames::geovals_testreference_logP})));
      profile.set<float>
        (addOPSPrefix(ufo::VariableNames::modellevels_ExnerP_rho_derived),
         std::move(profile.getGeoVaLVector
                   (oops::Variable{ufo::VariableNames::geovals_testreference_ExnerP_rho})));
      profile.set<float>
        (addOPSPrefix(ufo::VariableNames::modellevels_ExnerP_derived),
         std::move(profile.getGeoVaLVector
                   (oops::Variable{ufo::VariableNames::geovals_testreference_ExnerP})));
    } else {
      // Create a copy here because the vector will be used later in the routine.
      std::vector <float> logP = profile.get<float>(ufo::VariableNames::LogP_derived);
      profile.set<float>
        (addOPSPrefix(ufo::VariableNames::LogP_derived),
         std::move(logP));
      std::vector <float> bigPgaps = profile.get<float>(ufo::VariableNames::bigPgaps_derived);
      profile.set<float>
        (addOPSPrefix(ufo::VariableNames::bigPgaps_derived),
         std::move(bigPgaps));
    }
  }

}  // namespace ufo
