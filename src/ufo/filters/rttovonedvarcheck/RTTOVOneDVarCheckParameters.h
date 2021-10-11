/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_RTTOVONEDVARCHECK_RTTOVONEDVARCHECKPARAMETERS_H_
#define UFO_FILTERS_RTTOVONEDVARCHECK_RTTOVONEDVARCHECKPARAMETERS_H_

#include <string>
#include <vector>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/filters/FilterParametersBase.h"
#include "ufo/rttov/ObsRadianceRTTOVParameters.h"

namespace ufo {

/// Parameters controlling the operation of the ObsBoundsCheck filter.
class RTTOVOneDVarCheckParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(RTTOVOneDVarCheckParameters, FilterParametersBase)

 public:
  /// Path to the b-matrix file
  oops::RequiredParameter<std::string> BMatrix{"BMatrix", this};

  /// Path to the r-matrix file
  oops::RequiredParameter<std::string> RMatrix{"RMatrix", this};

  /// Specify the number of levels in the profiles
  oops::RequiredParameter<int> NLevels{"nlevels", this};

  /// List of the retrieval variables these will need to match the b-matrix file
  oops::RequiredParameter<std::vector<std::string>>
                                 RetrievalVariables{"retrieval variables", this};

  /// Options required for the forward model - RTTOV
  oops::RequiredParameter<ObsRadianceRTTOVParameters> ModOptions{"ModOptions", this};

  /// Specify the forward model to use - currently only RTTOV
  /// is available
  oops::Parameter<std::string> ModName{"ModName", "RTTOV", this};

  /// Is qtotal being used instead of separate q, clw, ciw
  oops::Parameter<bool> QTotal{"qtotal", false, this};

  /// Choose whether to split rain in qsplit routine
  oops::Parameter<bool> UseQtSplitRain{"UseQtSplitRain", true, this};

  /// Make sure profile is setup to use RTTOV-Scatt
  oops::Parameter<bool> RTTOVMWScattSwitch{"RTTOVMWScattSwitch", false, this};

  /// Use total ice option for RTTOV-Scatt
  oops::Parameter<bool> RTTOVUseTotalIce{"RTTOVUseTotalIce", true, this};

  /// Use the Marquardt-Levenberg minimizer - default is Newton
  oops::Parameter<bool> UseMLMinimization{"UseMLMinimization", false, this};

  /// Use cost function to determine when a profile has converged
  /// the default convergence options tests the absolute difference in the
  /// profile between iterations
  oops::Parameter<bool> UseJforConvergence{"UseJforConvergence", false, this};

  /// Use liquid water in the q saturation calculations
  oops::Parameter<bool> UseRHwaterForQC{"UseRHwaterForQC", true, this};

  /// Reset low level temperatures over seaice and cold, low land
  oops::Parameter<bool> UseColdSurfaceCheck{"UseColdSurfaceCheck", false, this};

  /// Output the LWP if the profile converges
  oops::Parameter<bool> Store1DVarLWP{"Store1DVarLWP", false, this};

  /// Turn on extra diagnostics
  oops::Parameter<bool> FullDiagnostics{"FullDiagnostics", false, this};

  /// Maximum number of iterations to perform in minimization
  oops::Parameter<int> Max1DVarIterations{"Max1DVarIterations", 7, this};

  /// Integer to select convergence option. 1= percentage change in cost tested between iterations
  /// otherwise = absolute change in cost tested between iterations
  oops::Parameter<int> JConvergenceOption{"JConvergenceOption", 1, this};

  /// Choose which iteration to start checking the liquid water path
  oops::Parameter<int> IterNumForLWPCheck{"IterNumForLWPCheck", 2, this};

  /// Maximum number of iterations for internal Marquardt-Levenberg loop
  oops::Parameter<int> MaxMLIterations{"MaxMLIterations", 7, this};

  /// Starting observation to run through 1d-var, subsetting for testing
  oops::Parameter<int> StartOb{"StartOb", 0, this};

  /// Final observation to run through 1d-var, subsetting for testing
  oops::Parameter<int> FinishOb{"FinishOb", 0, this};

  /// Check all the retrieved brightness temperatures are within a factor * error of the
  /// observed and bias corrected BTs.  If this value is less than 0.0 this check is
  /// not performed
  oops::Parameter<double> RetrievedErrorFactor{"RetrievedErrorFactor", 4.0, this};

  /// Convergence factor used when the absolute difference in the profile is used
  /// to determine convergence.
  oops::Parameter<double> ConvergenceFactor{"ConvergenceFactor", 0.4, this};

  /// Cost threshold for convergence check when cost function value is used for convergence
  oops::Parameter<double> CostConvergenceFactor{"CostConvergenceFactor", 0.01, this};

  /// Default emissivity value to use over land
  oops::Parameter<double> EmissLandDefault{"EmissLandDefault", 0.95, this};

  /// Default emissivity value to use over seaice
  oops::Parameter<double> EmissSeaIceDefault{"EmissSeaIceDefault", 0.92, this};

  /// Default eigen value path is blank but needs to be present if using PC emiss
  oops::Parameter<std::string> EmisEigVecPath{"EmisEigVecPath", "", this};

  /// Default emis atlas path is blank
  oops::Parameter<std::string> EmisAtlas{"EmisAtlas", "", this};
};

}  // namespace ufo

#endif  // UFO_FILTERS_RTTOVONEDVARCHECK_RTTOVONEDVARCHECKPARAMETERS_H_
