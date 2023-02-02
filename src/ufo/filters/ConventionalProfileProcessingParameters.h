/*
 * (C) Crown copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_CONVENTIONALPROFILEPROCESSINGPARAMETERS_H_
#define UFO_FILTERS_CONVENTIONALPROFILEPROCESSINGPARAMETERS_H_

#include <map>
#include <string>
#include <vector>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"

#include "ufo/filters/FilterParametersBase.h"

#include "ufo/profile/DataHandlerParameters.h"

#include "ufo/utils/Constants.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"
#include "ufo/utils/ProbabilityOfGrossErrorParameters.h"

namespace ufo {

  /// \brief Options controlling the operation of the ConventionalProfileProcessing filter.
  class ConventionalProfileProcessingParameters : public FilterParametersBase {
    OOPS_CONCRETE_PARAMETERS(ConventionalProfileProcessingParameters, FilterParametersBase)

   public:  // variables
    /// @name Generic parameters
    /// @{

    /// List of checks to perform
    oops::Parameter<std::vector<std::string>> Checks {"Checks", {}, this};

    /// Print station ID
    oops::Parameter<bool> PrintStationID {"PrintStationID", false, this};

    /// @}

    /// @name Standard level-related parameters
    /// @{

    /// Min P for finding standard levels (Pa)
    oops::Parameter<float> FS_MinP {"FS_MinP", 0.0, this};

    /// Standard Levels (hPa)
    oops::Parameter<std::vector<float>> StandardLevels{"StandardLevels",
        {1000, 925, 850, 700, 500, 400, 300, 250, 200, 150, 100, 70, 50, 30, 20, 10, 7, 3, 2, 1},
        this};

    /// @}

    /// @name Basic check parameters
    /// @{

    /// Skip basic checks. Should only be set to true in specific circumstances
    /// (e.g. for the relative humidity QC check;
    /// the OPS RH check routine does not apply the basic checks).
    oops::Parameter<bool> BChecks_Skip {"BChecks_Skip", false, this};

    /// Minimum value of pressure (Pa)
    oops::Parameter<float> BChecks_minValidP {"BChecks_minValidP", 0.0, this};

    /// Maximum value of pressure (Pa)
    oops::Parameter<float> BChecks_maxValidP {"BChecks_maxValidP", 110.0e3, this};

    /// Set flags for failed basic checks?
    oops::Parameter<bool> flagBasicChecksFail {"flagBasicChecksFail", true, this};

    /// @}

    /// @name Same P/different T check parameters
    /// @{

    /// Threshold used for same P/different T check (K)
    oops::Parameter<float> SPDTCheck_TThresh {"SPDTCheck_TThresh", 0.0, this};

    /// @}

    /// @name Sign check parameters
    /// @{

    /// Threshold used for Pstar difference in sign check (Pa)
    oops::Parameter<float> SCheck_PstarThresh {"SCheck_PstarThresh", 1000.0, this};

    /// Threshold used for |tObs - tBkg| in sign check (K)
    oops::Parameter<float> SCheck_tObstBkgThresh {"SCheck_tObstBkgThresh", 5.0, this};

    /// Tolerance used for sign check (K)
    oops::Parameter<float> SCheck_ProfileSignTol {"SCheck_ProfileSignTol", 100.0, this};

    /// P threshold over which to print large T differences (Pa)
    oops::Parameter<float> SCheck_PrintLargeTThresh {"SCheck_PrintLargeTThresh", 1000.0, this};

    /// Correct tObs in the sign check?
    oops::Parameter<bool> SCheck_CorrectT {"SCheck_CorrectT", true, this};

    /// @}

    /// @name Unstable layer check parameters
    /// @{

    /// Min P for unstable layer/superadiabat check (Pa)
    oops::Parameter<float> ULCheck_MinP {"ULCheck_MinP", 0.0, this};

    /// Bottom pressure threshold for unstable layer/superadiabat check (Pa)
    oops::Parameter<float> ULCheck_PBThresh {"ULCheck_PBThresh", 10000.0, this};

    /// Tolerance for unstable layer/superadiabat check (K)
    oops::Parameter<float> ULCheck_SuperadiabatTol {"ULCheck_SuperadiabatTol", -1.0, this};

    /// @}

    /// @name Interpolation check parameters
    /// @{

    /// Initial 'big gap' for interpolation check (hPa)
    oops::Parameter<float> ICheck_BigGapInit {"ICheck_BigGapInit", 1000.0, this};

    /// Pressure threshold for T tolerance relaxation
    oops::Parameter<float> ICheck_TolRelaxPThresh {"ICheck_TolRelaxPThresh", 50000.0, this};

    /// T tolerance relaxation factor
    oops::Parameter<float> ICheck_TolRelax {"ICheck_TolRelax", 1.0, this};

    /// Tolerance for interpolation check (K)
    oops::Parameter<float> ICheck_TInterpTol {"ICheck_TInterpTol", 1.0, this};

    /// Big gaps (hPa) used in interpolation check
    oops::Parameter<std::vector<float>> BigGaps{"ICheck_BigGaps",
        {500, 500, 500, 500, 100, 100, 100, 100,
            50, 50, 50, 50, 10, 10, 10, 10, 10, 10, 10, 10}, this};

    /// @}

    /// @name Hydrostatic check parameters

    /// @{

    /// Surface P threshold for hydrostatic check (Pa)
    oops::Parameter<float> HCheck_SurfacePThresh {"HCheck_SurfacePThresh", 10000.0, this};

    // A variety of thresholds used in the hydrostatic check
    oops::Parameter<float> HCheck_ETolMult {"HCheck_ETolMult", 0.5, this};
    oops::Parameter<float> HCheck_ETolMax {"HCheck_ETolMax", 1.0, this};
    oops::Parameter<float> HCheck_ETolMaxPThresh {"HCheck_ETolMaxPThresh", 50000.0, this};
    oops::Parameter<float> HCheck_ETolMaxLarger {"HCheck_ETolMaxLarger", 1.0, this};
    oops::Parameter<float> HCheck_ETolMin {"HCheck_ETolMin", 1.0, this};
    oops::Parameter<float> HCheck_EThresh {"HCheck_EThresh", 100.0, this};
    oops::Parameter<float> HCheck_EThreshB {"HCheck_EThreshB", 100.0, this};
    oops::Parameter<float> HCheck_ESumThresh {"HCheck_ESumThresh", 50.0, this};
    oops::Parameter<float> HCheck_MinAbsEThresh {"HCheck_MinAbsEThresh", 10.0, this};
    oops::Parameter<float> HCheck_ESumThreshLarger {"HCheck_ESumThreshLarger", 100.0, this};
    oops::Parameter<float> HCheck_MinAbsEThreshLarger {"HCheck_MinAbsEThreshLarger", 100.0, this};
    oops::Parameter<float> HCheck_CorrThresh {"HCheck_CorrThresh", 5.0, this};
    oops::Parameter<float> HCheck_ESumNextThresh {"HCheck_ESumNextThresh", 50.0, this};
    oops::Parameter<float> HCheck_MinAbsEThreshT {"HCheck_MinAbsEThreshT", 10.0, this};
    oops::Parameter<float> HCheck_CorrDiffThresh {"HCheck_CorrDiffThresh", 10.0, this};
    oops::Parameter<float> HCheck_CorrMinThresh {"HCheck_CorrMinThresh", 1.0, this};

    /// Correct zObs in the hydrostatic check?
    oops::Parameter<bool> HCheck_CorrectZ {"HCheck_CorrectZ", true, this};

    /// Hydrostatic error descriptions
    oops::Parameter<std::vector<std::string>> HydDesc{"HydDesc",
        {"Hyd: OK", "Hyd: Z err", "Hyd: T err",
            "Hyd: T/Z err", "Hyd: T/Z Bot",
            "Hyd: T/Z Top", "Hyd: Z upward", "Hyd: ?????"},
        this};

    /// @}

    /// @name Wind speed interpolation check parameters
    /// @{

    /// Squared tolerance for identical pressure in wind speed interpolation check (m^2 s^-2)
    oops::Parameter<float> UICheck_TInterpIdenticalPTolSq
      {"UICheck_TInterpIdenticalPTolSq", 0.0, this};

    /// Squared tolerance for wind speed interpolation check (m^2 s^-2)
    oops::Parameter<float> UICheck_TInterpTolSq {"UICheck_TInterpTolSq", 0.0, this};

    /// Big gap (Pa) used at lowest pressures in wind speed interpolation check
    oops::Parameter<float> UICheck_BigGapLowP {"UICheck_BigGapLowP", 500.0, this};

    /// Big gaps (Pa) used in wind speed interpolation check.
    /// This vector must be the same length as UICheck_BigGapsPThresh.
    oops::Parameter<std::vector<float>> UICheck_BigGaps{"UICheck_BigGaps",
        {50000.0, 10000.0, 5000.0, 1000.0}, this};

    /// Big gap thresholds (Pa) used in wind speed interpolation check.
    /// This vector must be the same length as UICheck_BigGaps.
    oops::Parameter<std::vector<float>> UICheck_BigGapsPThresh{"UICheck_BigGapsPThresh",
        {100000.0, 50000.0, 10000.0, 5000.0}, this};

    /// @}

    /// @name RH check parameters
    /// @{

    /// Initial value of minimum temperature (K)
    oops::Parameter<float> RHCheck_TminInit {"RHCheck_TminInit", 400.0, this};

    /// Tolerance for high level check of relative humidity (%)
    oops::Parameter<float> RHCheck_SondeRHHiTol {"RHCheck_SondeRHHiTol", 0.0, this};

    /// Threshold for pressure when setting up arrays (Pa)
    oops::Parameter<float> RHCheck_PressInitThresh {"RHCheck_PressInitThresh", 500.0, this};

    /// Threshold for pressure (Pa)
    oops::Parameter<float> RHCheck_PressThresh {"RHCheck_PressThresh", 500.0, this};

    /// Threshold for pressure difference relative to level 0 (Pa)
    oops::Parameter<float> RHCheck_PressDiff0Thresh {"RHCheck_PressDiff0Thresh", 50.0, this};

    /// Threshold for dew point temperature difference (K)
    oops::Parameter<float> RHCheck_tdDiffThresh {"RHCheck_tdDiffThresh", 5.0, this};

    /// Threshold for relative humidity (%)
    oops::Parameter<float> RHCheck_RHThresh {"RHCheck_RHThresh", 75.0, this};

    /// Threshold for pressure difference between adjacent levels (Pa)
    oops::Parameter<float> RHCheck_PressDiffAdjThresh {"RHCheck_PressDiffAdjThresh", 50.0, this};

    /// Threshold for minimum relative humidity (%)
    oops::Parameter<float> RHCheck_MinRHThresh {"RHCheck_MinRHThresh", 75.0, this};

    /// Upper threshold for Tmin in moisture check
    oops::Parameter<float> RHCheck_TminThresh {"RHCheck_TminThresh", 200.0, this};

    /// Lower threshold for temperature in moisture check
    oops::Parameter<float> RHCheck_TempThresh {"RHCheck_TempThresh", 250.0, this};

    /// @}

    /// @name Time check parameters
    /// @{

    /// Threshold relative to surface pressure for rejecting levels (hPa)
    oops::Parameter<float> TimeCheck_SondeLaunchWindRej {"TimeCheck_SondeLaunchWindRej", 0.0, this};

    /// @}

    /// @name Background check (T, RH, UV) parameters
    /// @{

    /// Prior probability of 'bad' observations for T
    oops::Parameter<float> BkCheck_PdBad_t {"BkCheck_PdBad_t", 0.05, this};

    /// Prior probability of 'bad' observations for RH
    oops::Parameter<float> BkCheck_PdBad_rh {"BkCheck_PdBad_rh", 0.05, this};

    /// Prior probability of 'bad' observations for u and v
    oops::Parameter<float> BkCheck_PdBad_uv {"BkCheck_PdBad_uv", 0.001, this};

    /// Observations with a latitude smaller than this value (both N and S)
    /// are taken to be in the tropics.
    oops::Parameter<float> BkCheck_Psplit_latitude_tropics
      {"BkCheck_Psplit_latitude_tropics", 30.0, this};

    /// Pressure threshold above which extra representivity error occurs in extratropics (Pa).
    oops::Parameter<float> BkCheck_Psplit_extratropics
      {"BkCheck_Psplit_extratropics", 50000.0, this};

    /// Pressure threshold above which extra representivity error occurs in tropics (Pa).
    oops::Parameter<float> BkCheck_Psplit_tropics {"BkCheck_Psplit_tropics", 10000.0, this};

    /// Error inflation factor below Psplit
    oops::Parameter<float> BkCheck_ErrorInflationBelowPsplit
      {"BkCheck_ErrorInflationBelowPsplit", 1.0, this};

    /// Error inflation factor above Psplit
    oops::Parameter<float> BkCheck_ErrorInflationAbovePsplit
      {"BkCheck_ErrorInflationAbovePsplit", 1.0, this};

    /// Maximum error variance for RH
    oops::Parameter<float> BkCheck_ErrVarMax_rh {"BkCheck_ErrVarMax_rh", 500.0, this};

    /// Pressure thresholds for setting z background errors and 'bad' observation PGE.
    /// This vector must be the same length as BkCheck_zBkgErrs and BkCheck_zBadPGEs.
    oops::Parameter<std::vector<float>> BkCheck_PlevelThresholds {"BkCheck_PlevelThresholds",
        {1000.0, 500.0, 100.0, 50.0, 10.0, 5.0, 1.0, 0.0},
        this};

    /// List of z background errors that are assigned based on pressure.
    /// This vector must be the same length as BkCheck_PlevelThresholds and BkCheck_zBadPGEs.
    oops::Parameter<std::vector<float>> BkCheck_zBkgErrs {"BkCheck_zBkgErrs",
        {10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0},
        this};

    /// List of z PGEs for 'bad' observations that are assigned based on pressure.
    /// This vector must be the same length as BkCheck_PlevelThresholds and BkCheck_zBkgErrs.
    oops::Parameter<std::vector<float>> BkCheck_zBadPGEs {"BkCheck_zBadPGEs",
        {0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01},
        this};

    /// @}

    /// @name Profile averaging parameters
    /// @{

    /// Factor used to determine big gaps for sondes
    /// (dimensionless; multiplied by log(10)).
    oops::Parameter<float> AvgP_SondeGapFactor {"AvgP_SondeGapFactor", 1.0, this};

    /// Factor used to determine big gaps for wind profilers
    /// (dimensionless; multiplied by log(10)).
    oops::Parameter<float> AvgP_WinProGapFactor {"AvgP_WinProGapFactor", 1.0, this};

    /// Minimum value of denominator used when computing big gaps
    /// (dimensionless; equal to log (pressure threshold / hPa)).
    oops::Parameter<float> AvgP_GapLogPDiffMin {"AvgP_GapLogPDiffMin", std::log(5.0), this};

    /// Minimum fraction of a model layer that must have been covered (in the vertical coordinate)
    /// by observed values in order for temperature to be averaged onto that layer.
    oops::Parameter<float> AvgT_SondeDZFraction {"AvgT_SondeDZFraction", 0.5, this};

    /// Probability of gross error threshold above which rejection flags are set
    /// in the temperature averaging routine.
    oops::Parameter<float> AvgT_PGEskip {"AvgT_PGEskip", 0.9, this};

    /// Minimum fraction of a model layer that must have been covered (in the vertical coordinate)
    /// by observed values in order for wind speed to be averaged onto that layer.
    oops::Parameter<float> AvgU_SondeDZFraction {"AvgU_SondeDZFraction", 0.5, this};

    /// Probability of gross error threshold above which rejection flags are set
    /// in the wind speed averaging routine.
    oops::Parameter<float> AvgU_PGEskip {"AvgU_PGEskip", 0.9, this};

    /// Probability of gross error threshold above which rejection flags are set
    /// in the relative humidity averaging routine.
    oops::Parameter<float> AvgRH_PGEskip {"AvgRH_PGEskip", 0.9, this};

    /// Minimum fraction of a model layer that must have been covered (in the vertical coordinate)
    /// by observed values in order for relative humidity to be averaged onto that layer.
    oops::Parameter<float> AvgRH_SondeDZFraction {"AvgRH_SondeDZFraction", 0.5, this};

    /// Perform interpolation or averaging of relative humidity observations?
    oops::Parameter<bool> AvgRH_Interp {"AvgRH_Interp", true, this};

    /// Default average temperature threshold below which average relative humidity
    /// observations are rejected (degrees C).
    oops::Parameter<float> AvgRH_AvgTThreshold {"AvgRH_AvgTThreshold", -40.0, this};

    /// Custom average temperature thresholds below which average relative humidity
    /// observations are rejected (degrees C).
    /// These thresholds are stored in a map with keys equal to the WMO codes for
    /// radiosonde instrument types and values equal to the custom thresholds.
    ///
    /// The full list of codes can be found in "WMO Manual on Codes -
    /// International Codes, Volume I.2, Annex II to the WMO Technical Regulations:
    /// Part C - Common Features to Binary and Alphanumeric Codes"
    /// (available at https://library.wmo.int/?lvl=notice_display&id=10684).
    ///
    /// The default list of custom thresholds applies to the following codes:
    /// - Types 37, 52, 60-63, 66-67 = Vaisala RS80,
    /// - Types 71-74, 78 = Vaisala RS90,
    /// - Types 79-81 = Vaisala RS92.
    ///
    /// To customise this list in the yaml file, please note the following information from
    /// oops/util/parameter/ParameterTraits.h:
    /// Owing to a bug in the eckit YAML parser, maps need to be written in the JSON style,
    /// with keys quoted. Example:
    ///   my_int_to_float_map: {"1": 123, "2": 321}
    oops::Parameter<std::map<int, float>> AvgRH_InstrTThresholds
    {"AvgRH_InstrTThresholds",
        {{37, -60.0}, {52, -60.0}, {60, -60.0}, {61, -60.0},
         {62, -60.0}, {63, -60.0}, {66, -60.0}, {67, -60.0},
         {71, -80.0}, {72, -80.0}, {73, -80.0}, {74, -80.0},
         {78, -80.0}, {79, -80.0}, {80, -80.0}, {81, -80.0}
        }, this};

    /// @}

    /// @name Background error names
    /// @{

    oops::Parameter<std::string> bkgErrGroup {"background error group", "ObsDiag", this};
    oops::Parameter<std::string> bkgErrName_air_temperature
      {"T background error name", "air_temperature_background_error", this};
    oops::Parameter<std::string> bkgErrName_eastward_wind
      {"u background error name", "eastward_wind_background_error", this};
    oops::Parameter<std::string> bkgErrName_northward_wind
      {"v background error name", "northward_wind_background_error", this};
    oops::Parameter<std::string> bkgErrName_relative_humidity
      {"RH background error name", "relative_humidity_background_error", this};

    /// @}

    /// @name OPS comparison parameters
    /// @{

    /// Compare with OPS values?
    oops::Parameter<bool> compareWithOPS {"compareWithOPS", false, this};

    /// Tolerance for absolute difference comparisions
    oops::Parameter<float> Comparison_Tol {"Comparison_Tol", 0.1, this};

    /// @}

    /// @name Parameters classes
    /// @{

    /// Parameters related to profile data handler
    DataHandlerParameters DHParameters{this};

    /// Parameters related to PGE calculations
    ProbabilityOfGrossErrorParameters PGEParameters{this};

    /// @}
  };
}  // namespace ufo

#endif  // UFO_FILTERS_CONVENTIONALPROFILEPROCESSINGPARAMETERS_H_

