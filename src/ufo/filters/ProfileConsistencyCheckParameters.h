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

#include "ufo/utils/Constants.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace eckit {
  class Configuration;
}

namespace ufo {

  /// \brief Options controlling the operation of the ProfileConsistencyChecks filter.
  class ProfileConsistencyCheckParameters : public oops::Parameters {
   public:  // functions
    /// Determine whether a variable group is optional or not.
    bool getOptional(const std::string &groupname) const
    {
      bool optional = false;
      if (std::find(groups_optional.value().begin(), groups_optional.value().end(), groupname)
          != groups_optional.value().end())
        optional = true;
      return optional;
    }

    /// Determine number of entries per profile for a variable group.
    size_t getEntriesPerProfile(const std::string &groupname) const
    {
      size_t entriesPerProfile = -1;
      if (std::find(groups_singlevalue.value().begin(), groups_singlevalue.value().end(), groupname)
          != groups_singlevalue.value().end())
        entriesPerProfile = 1;
      if (std::find(groups_modellevels.value().begin(), groups_modellevels.value().end(), groupname)
          != groups_modellevels.value().end())
        entriesPerProfile = static_cast<size_t> (num_modellevels);
      return entriesPerProfile;
    }

   public:  // variables
    //=== Generic parameters ===//

    /// Maximum number of profile levels to be processed (a legacy of the OPS code).
    /// No maximum is assigned if this parameter is not specified.
    oops::OptionalParameter<int> maxlev {"maxlev", this};

    /// List of checks to perform
    oops::Parameter<std::vector<std::string>> Checks {"Checks", {}, this};

    /// If not sorting observations, ensure number of profiles is consistent
    oops::Parameter<bool> ValidateTotalNumProf {"ValidateTotalNumProf", true, this};

    /// Print station ID
    oops::Parameter<bool> PrintStationID {"PrintStationID", false, this};

    //=== Parameters relating to the combination of multiple check results ===//

    /// Number of errors that cause the observation to have failed
    oops::Parameter<int> nErrorsFail {"nErrorsFail", 8, this};

    ///=== Standard level-related parameters ===//

    /// Min P for finding standard levels (Pa)
    oops::Parameter<float> FS_MinP {"FS_MinP", 0.0, this};

    /// Standard Levels (hPa)
    oops::Parameter<std::vector<float>> StandardLevels{"StandardLevels",
        {1000, 925, 850, 700, 500, 400, 300, 250, 200, 150, 100, 70, 50, 30, 20, 10, 7, 3, 2, 1},
        this};

    //=== Basic check parameters ===//

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

    //=== Same P/different T check parameters ===//

    /// Threshold used for same P/different T check (K)
    oops::Parameter<float> SPDTCheck_TThresh {"SPDTCheck_TThresh", 1.0, this};

    //=== Sign check parameters ===//

    /// Threshold used for Pstar difference in sign check (Pa)
    oops::Parameter<float> SCheck_PstarThresh {"SCheck_PstarThresh", 5000.0, this};

    /// Threshold used for |tObs - tBkg| in sign check (K)
    oops::Parameter<float> SCheck_tObstBkgThresh {"SCheck_tObstBkgThresh", 20.0, this};

    /// Tolerance used for sign check (K)
    oops::Parameter<float> SCheck_ProfileSignTol {"SCheck_ProfileSignTol", 5.0, this};

    /// P threshold over which to print large T differences (Pa)
    oops::Parameter<float> SCheck_PrintLargeTThresh {"SCheck_PrintLargeTThresh", 1000.0, this};

    /// Correct tObs in the sign check?
    oops::Parameter<bool> SCheck_CorrectT {"SCheck_CorrectT", true, this};

    //=== Unstable layer check parameters ===//

    /// Min P for unstable layer/superadiabat check (Pa)
    oops::Parameter<float> ULCheck_MinP {"ULCheck_MinP", 0.0, this};

    /// Bottom pressure threshold for unstable layer/superadiabat check (Pa)
    oops::Parameter<float> ULCheck_PBThresh {"ULCheck_PBThresh", 5000.0, this};

    /// Tolerance for unstable layer/superadiabat check (K)
    oops::Parameter<float> ULCheck_SuperadiabatTol {"ULCheck_SuperadiabatTol", -2.0, this};

    //=== Interpolation check parameters ===//

    /// Initial 'big gap' for interpolation check (hPa)
    oops::Parameter<float> ICheck_BigGapInit {"ICheck_BigGapInit", 1000.0, this};

    /// Pressure threshold for T tolerance relaxation
    oops::Parameter<float> ICheck_TolRelaxPThresh {"ICheck_TolRelaxPThresh", 30000.0, this};

    /// T tolerance relaxation factor
    oops::Parameter<float> ICheck_TolRelax {"ICheck_TolRelax", 1.5, this};

    /// Tolerance for interpolation check (K)
    oops::Parameter<float> ICheck_TInterpTol {"ICheck_TInterpTol", 2.0, this};

    /// Big gaps (hPa) used in interpolation check
    oops::Parameter<std::vector<float>> BigGaps{"ICheck_BigGaps",
        {150, 150, 150, 150, 100, 100, 100, 75,
            75, 50, 50, 20, 20, 20, 10, 10, 10, 10, 10, 10}};

    //=== Hydrostatic check parameters ===//

    /// Surface P threshold for hydrostatic check (Pa)
    oops::Parameter<float> HCheck_SurfacePThresh {"HCheck_SurfacePThresh", 15100.0, this};

    // A variety of thresholds used in the hydrostatic check
    oops::Parameter<float> HCheck_ETolMult {"HCheck_ETolMult", 0.375, this};
    oops::Parameter<float> HCheck_ETolMax {"HCheck_ETolMax", 50.0, this};
    oops::Parameter<float> HCheck_ETolMaxPThresh {"HCheck_ETolMaxPThresh", 40100.0, this};
    oops::Parameter<float> HCheck_ETolMaxLarger {"HCheck_ETolMaxLarger", 80.0, this};
    oops::Parameter<float> HCheck_ETolMin {"HCheck_ETolMin", 30.0, this};
    oops::Parameter<float> HCheck_EThresh {"HCheck_EThresh", 15.0, this};
    oops::Parameter<float> HCheck_EThreshB {"HCheck_EThreshB", 15.0, this};
    oops::Parameter<float> HCheck_ESumThresh {"HCheck_ESumThresh", 30.0, this};
    oops::Parameter<float> HCheck_MinAbsEThresh {"HCheck_MinAbsEThresh", 20.0, this};
    oops::Parameter<float> HCheck_ESumThreshLarger {"HCheck_ESumThreshLarger", 60.0, this};
    oops::Parameter<float> HCheck_MinAbsEThreshLarger {"HCheck_MinAbsEThreshLarger", 200.0, this};
    oops::Parameter<float> HCheck_CorrThresh {"HCheck_CorrThresh", 10.0, this};
    oops::Parameter<float> HCheck_ESumNextThresh {"HCheck_ESumNextThresh", 30.0, this};
    oops::Parameter<float> HCheck_MinAbsEThreshT {"HCheck_MinAbsEThreshT", 15.0, this};
    oops::Parameter<float> HCheck_CorrDiffThresh {"HCheck_CorrDiffThresh", 5.0, this};
    oops::Parameter<float> HCheck_CorrMinThresh {"HCheck_CorrMinThresh", 4.0, this};

    /// Correct zObs in the hydrostatic check?
    oops::Parameter<bool> HCheck_CorrectZ {"HCheck_CorrectZ", true, this};

    /// Hydrostatic error descriptions
    oops::Parameter<std::vector<std::string>> HydDesc{"HydDesc",
        {"Hyd: OK", "Hyd: Z err", "Hyd: T err",
            "Hyd: T/Z err", "Hyd: T/Z Bot",
            "Hyd: T/Z Top", "Hyd: Z upward", "Hyd: ?????"},
        this};

    //=== Wind speed interpolation check parameters ===//

    /// Squared tolerance for identical pressure in wind speed interpolation check (m^2 s^-2)
    oops::Parameter<float> UICheck_TInterpIdenticalPTolSq
      {"UICheck_TInterpIdenticalPTolSq", 4.0, this};

    /// Squared tolerance for wind speed interpolation check (m^2 s^-2)
    oops::Parameter<float> UICheck_TInterpTolSq {"UICheck_TInterpTolSq", 64.0, this};

    /// Big gap (Pa) used at lowest pressures in wind speed interpolation check
    oops::Parameter<float> UICheck_BigGapLowP {"UICheck_BigGapLowP", 1000.0, this};

    /// Big gaps (Pa) used in wind speed interpolation check
    oops::Parameter<std::vector<float>> UICheck_BigGaps{"UICheck_BigGaps",
        {15000.0, 10000.0, 7500.0, 5000.0, 2000.0}};

    /// Big gap thresholds (Pa) used in wind speed interpolation check
    oops::Parameter<std::vector<float>> UICheck_BigGapsPThresh{"UICheck_BigGapsPThresh",
        {65000.0, 27500.0, 17500.0, 8500.0, 2500.0}};

    //=== RH check parameters ===//

    /// Initial value of minimum temperature (K)
    oops::Parameter<float> RHCheck_TminInit {"RHCheck_TminInit", 400.0, this};

    /// Tolerance for high level check of relative humidity (%)
    oops::Parameter<float> RHCheck_SondeRHHiTol {"RHCheck_SondeRHHiTol", 20.0, this};

    /// Threshold for pressure when setting up arrays (Pa)
    oops::Parameter<float> RHCheck_PressInitThresh {"RHCheck_PressInitThresh", 100.0, this};

    /// Threshold for pressure (Pa)
    oops::Parameter<float> RHCheck_PressThresh {"RHCheck_PressThresh", 400.0, this};

    /// Threshold for pressure difference relative to level 0 (Pa)
    oops::Parameter<float> RHCheck_PressDiff0Thresh {"RHCheck_PressDiff0Thresh", 100.0, this};

    /// Threshold for dew point temperature difference (K)
    oops::Parameter<float> RHCheck_tdDiffThresh {"RHCheck_tdDiffThresh", 2.0, this};

    /// Threshold for relative humidity (%)
    oops::Parameter<float> RHCheck_RHThresh {"RHCheck_RHThresh", 90.0, this};

    /// Threshold for pressure difference between adjacent levels (Pa)
    oops::Parameter<float> RHCheck_PressDiffAdjThresh {"RHCheck_PressDiffAdjThresh", 20.0, this};

    /// Threshold for minimum relative humidity (%)
    oops::Parameter<float> RHCheck_MinRHThresh {"RHCheck_MinRHThresh", 85.0, this};

    /// Upper threshold for Tmin in moisture check
    oops::Parameter<float> RHCheck_TminThresh {"RHCheck_TminThresh", 223.15, this};

    /// Lower threshold for temperature in moisture check
    oops::Parameter<float> RHCheck_TempThresh {"RHCheck_TempThresh", 233.15, this};

    //=== OPS comparison parameters ===//

    /// Tolerance for absolute difference comparisions
    oops::Parameter<float> Comparison_Tol {"Comparison_Tol", 0.1, this};

    /// Compare with OPS values?
    oops::Parameter<bool> compareWithOPS {"compareWithOPS", false, this};

    /// Groups of variables whose presence in the input sample is optional
    /// (if not present, all variables are initially set to zero)
    oops::Parameter<std::vector<std::string>> groups_optional
      {"groups_optional", {"Corrections", "Counters"}, this};

    /// Groups of variables which have one value per profile
    oops::Parameter<std::vector<std::string>> groups_singlevalue
      {"groups_singlevalue", {"Counters"}, this};

    /// Groups of variables with values on model levels
    oops::Parameter<std::vector<std::string>> groups_modellevels
      {"groups_modellevels", {}, this};

    /// Number of model levels
    oops::Parameter<int> num_modellevels {"num_modellevels", 70, this};
  };
}  // namespace ufo

#endif  // UFO_FILTERS_PROFILECONSISTENCYCHECKPARAMETERS_H_

