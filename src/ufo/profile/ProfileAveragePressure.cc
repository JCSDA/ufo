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

#include "ufo/utils/metoffice/MetOfficeObservationIDs.h"

namespace ufo {

  static ProfileCheckMaker<ProfileAveragePressure>
  makerProfileAveragePressure_("AveragePressure");

  ProfileAveragePressure::ProfileAveragePressure
  (const ProfileConsistencyCheckParameters &options)
    : ProfileCheckBase(options)
  {}

  void ProfileAveragePressure::logPressure(const std::vector <float> &pressures,
                                           std::vector <float> &logP)
  {
    logP = pressures;
    std::transform(logP.begin(), logP.end(), logP.begin(),
                   [this](float PLev){return PLev > 0 ? log(PLev) : missingValueFloat;});
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
    const float GapSize = log(10.0f);
    const float GapFactor = ObsType == ufo::MetOfficeObsIDs::AtmosphericProfile::WindProf ?
      options_.AvgP_WinProGapFactor.value() * GapSize :
      options_.AvgP_SondeGapFactor.value() * GapSize;
    const float GapLogPDiffMin = options_.AvgP_GapLogPDiffMin.value();
    // Subtracting log(100.0) from PLev effectively converts the pressure to hPa.
    std::transform(bigPgaps.begin(), bigPgaps.end(), bigPgaps.begin(),
                   [this, GapFactor, GapLogPDiffMin](float PLev){return PLev > 0.0 ?
                       GapFactor / std::max(log(PLev) - log(100.0f), GapLogPDiffMin) :
                       missingValueFloat;});
  }

  void ProfileAveragePressure::runCheck(ProfileDataHandler &profileDataHandler)
  {
    oops::Log::debug() << " Pressure transformation" << std::endl;

    const size_t numProfileLevels = profileDataHandler.getNumProfileLevels();
    const std::vector <float> &pressures =
      profileDataHandler.get<float>(ufo::VariableNames::obs_air_pressure);
    const std::vector <int> &ObsType =
      profileDataHandler.get<int>(ufo::VariableNames::ObsType);
    const std::vector <int> &ReportFlags =
      profileDataHandler.get<int>(ufo::VariableNames::qcflags_observation_report);

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
    const size_t numModelLevels_rho = options_.DHParameters.ModParameters.numModelLevels_rho();
    std::vector <float> geovals_logP_rho(numModelLevels_rho, missingValueFloat);  // log(pressure)
    std::vector <float> geovals_ExnerP_rho
      (numModelLevels_rho, missingValueFloat);  // Exner pressure

    // Initialise vectors on model theta levels.
    const size_t numModelLevels = options_.DHParameters.ModParameters.numModelLevels();
    std::vector <float> geovals_logP(numModelLevels, missingValueFloat);  // log(pressure)
    std::vector <float> geovals_ExnerP(numModelLevels, missingValueFloat);  // Exner pressure

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
        profileDataHandler.getGeoVaLVector(ufo::VariableNames::geovals_pressure_rho);
      logPressure(geovals_pressure_rho, geovals_logP_rho);
      ExnerPressure(geovals_pressure_rho, geovals_ExnerP_rho);

      // Calculate log(P) and Exner pressure for GeoVaLs on theta levels.
      const std::vector <float> &geovals_pressure =
        profileDataHandler.getGeoVaLVector(ufo::VariableNames::geovals_pressure);
      // Require these GeoVaLs to be present (unlike the case on rho levels).
      if (geovals_pressure.empty())
        throw eckit::BadValue("geovals_pressure is empty", Here());
      logPressure(geovals_pressure, geovals_logP);
      ExnerPressure(geovals_pressure, geovals_ExnerP);
    }

    // Store the transformed pressures on reported levels.
    profileDataHandler.set<float>(ufo::VariableNames::LogP_derived, std::move(logP));
    profileDataHandler.set<float>(ufo::VariableNames::bigPgaps_derived, std::move(bigPgaps));

    // Store the transformed pressures on model levels.
    profileDataHandler.set<float>(ufo::VariableNames::geovals_logP_rho_derived,
                                  std::move(geovals_logP_rho));
    profileDataHandler.set<float>(ufo::VariableNames::geovals_logP_derived,
                                  std::move(geovals_logP));
    profileDataHandler.set<float>(ufo::VariableNames::geovals_ExnerP_rho_derived,
                                  std::move(geovals_ExnerP_rho));
    profileDataHandler.set<float>(ufo::VariableNames::geovals_ExnerP_derived,
                                  std::move(geovals_ExnerP));
  }

  void ProfileAveragePressure::fillValidationData(ProfileDataHandler &profileDataHandler)
  {
    // Retrieve, then save, the OPS versions of the transformed pressures on model levels.
    profileDataHandler.set<float>
      ("OPS_" + std::string(ufo::VariableNames::geovals_logP_rho_derived),
       std::move(profileDataHandler.getGeoVaLVector
                 (ufo::VariableNames::geovals_logP_rho)));
    profileDataHandler.set<float>
      ("OPS_" + std::string(ufo::VariableNames::geovals_logP_derived),
       std::move(profileDataHandler.getGeoVaLVector
                 (ufo::VariableNames::geovals_logP)));
    profileDataHandler.set<float>
      ("OPS_" + std::string(ufo::VariableNames::geovals_ExnerP_rho_derived),
       std::move(profileDataHandler.getGeoVaLVector
                 (ufo::VariableNames::geovals_ExnerP_rho)));
    profileDataHandler.set<float>
      ("OPS_" + std::string(ufo::VariableNames::geovals_ExnerP_derived),
       std::move(profileDataHandler.getGeoVaLVector
                 (ufo::VariableNames::geovals_ExnerP)));
  }
}  // namespace ufo
