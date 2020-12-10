/*
 * (C) Crown copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_PROFILECONSISTENCYCHECKPARAMETERS_H_
#define UFO_FILTERS_PROFILECONSISTENCYCHECKPARAMETERS_H_

#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"

#include "ufo/profile/DataHandlerParameters.h"

#include "ufo/utils/Constants.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"
#include "ufo/utils/ProbabilityOfGrossErrorParameters.h"

namespace eckit {
  class Configuration;
}

namespace ufo {

  /// \brief Options controlling the operation of the ProfileConsistencyChecks filter.
  class ProfileConsistencyCheckParameters : public oops::Parameters {
    OOPS_CONCRETE_PARAMETERS(ProfileConsistencyCheckParameters, Parameters)

   public:  // variables
    /// @name Generic parameters
    /// @{

    /// List of checks to perform
    oops::Parameter<std::vector<std::string>> Checks {"Checks", {}, this};

    /// Print station ID
    oops::Parameter<bool> PrintStationID {"PrintStationID", false, this};

    /// Have the observation and model values been averaged onto model levels?
    oops::Parameter<bool> modellevels {"ModelLevels", false, this};

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

#endif  // UFO_FILTERS_PROFILECONSISTENCYCHECKPARAMETERS_H_

