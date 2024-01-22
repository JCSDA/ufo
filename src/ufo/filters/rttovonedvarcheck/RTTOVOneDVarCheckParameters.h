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
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/filters/FilterParametersBase.h"
#include "ufo/operators/rttov/ObsRadianceRTTOVParameters.h"

namespace ufo {

/// Parameters class for the surface emissivity variables
class SurfaceEmissivityParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SurfaceEmissivityParameters, Parameters)

 public:
  /// How to initialise the surface emissivity
  /// rttovtocalculate - rttov will calculate over all surfaces
  /// fixed - values of EmissSeaDefault, EmissLandDefault, EmissSeaIceDefault
  /// readfromdb - read values from the db that have previously been read in
  ///              and set calc_emiss to true where zero.
  /// readfromdbwitherror - read values from the db that have previously
  ///              been read in and set calc_emiss to true where zero.
  ///              Read in the emissivity error as well.
  /// principalcomponent - net ready yet
  oops::Parameter<std::string> type{"type", "fixed", this};

  /// Default emissivity value to use over land - zero means RTTOV will calculate
  /// used for type = fixed
  oops::Parameter<double> EmissSeaDefault{"EmissSeaDefault", 0.00, this};

  /// Default emissivity value to use over land - zero means RTTOV will calculate
  /// used for type = fixed
  oops::Parameter<double> EmissLandDefault{"EmissLandDefault", 0.95, this};

  /// Default emissivity value to use over seaice - zero means RTTOV will calculate
  /// used for type = fixed
  oops::Parameter<double> EmissSeaIceDefault{"EmissSeaIceDefault", 0.92, this};

  /// Location of emissivity values to be read from the database e.g.
  /// default is DerivedObsValue/emissivity_<chan>
  oops::Parameter<std::string> groupInObsSpace{"group in obs space", "DerivedObsValue", this};

  /// Default eigen value path is blank but needs to be present if using PC emiss
  /// not currently used to be implemented
  /// used for type = hyperspectralpc
  oops::Parameter<std::string> EmisEigVecPath{"EmisEigVecPath", "", this};

  /// Default emis atlas path is blank - not currently used to be implemented
  /// used for type = hyperspectralpc
  oops::Parameter<std::string> EmisAtlas{"EmisAtlas", "", this};

  /// Flag to decide if mwemiss retrieval needed
  oops::Parameter<bool> mwEmissRetrieval{"retrieve mw emissivity", false, this};

  /// Number of surface emissivity elements to be retrieved.  This will be checked
  /// against the number in the b-matrix.
  /// used when surface_emissivity is retrieved unless pcemiss is specified
  oops::Parameter<int> NumEmissElements{
      "number of surface emissivity retrieval elements", 5, this};

  /// Maps the correct emissivity element to the correct instrument channel.
  /// Must be of size NumEmissElements.
  /// used when surface_emissivity is retrieved unless pcemiss is specified
  oops::Parameter<std::vector<int>> EmissToChannelMap{
      "emissivity to channel mapping", {1, 2, 3, 16, 17}, this};

  /// Maps the instrument channels to the correct emissivity element used in the retrieval.
  /// used when surface_emissivity is retrieved unless pcemiss is specified
  oops::Parameter<std::vector<int>> ChannelToEmissMap{
      "channel to emissivity mapping", {1, 2, 3, 3, 3, 3, 3, 3, 3, 3,
                                        3, 3, 3, 3, 4, 4, 5, 5, 5, 5}, this};
};

/// Parameters controlling the operation of the RTTOVOneDVarCheck filter.
class RTTOVOneDVarCheckParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(RTTOVOneDVarCheckParameters, FilterParametersBase)

 public:
  /// Path to the b-matrix file
  oops::RequiredParameter<std::string> BMatrix{"BMatrix", this};

  /// Path to the r-matrix file
  oops::RequiredParameter<std::string> RMatrix{"RMatrix", this};

  /// Specify the number of levels in the profiles
  oops::RequiredParameter<int> NLevels{"nlevels", this};

  /// List of the retrieval variables that are in the GeoVaLs
  /// these will need to match the b-matrix file
  oops::RequiredParameter<std::vector<std::string>> RetrievalVariablesGeoVaLs{
      "retrieval variables from geovals", this};

  /// List of the retrieval variables that are not in the GeoVaLs
  /// these will need to match the b-matrix file
  oops::Parameter<std::vector<std::string>> RetrievalVariablesNotGeoVaLs{
      "retrieval variables not from geovals", {}, this};

  /// Options required for the forward model - RTTOV
  oops::RequiredParameter<ObsRadianceRTTOVParameters> ModOptions{"ModOptions", this};

  /// Specify the forward model to use - currently only RTTOV
  /// is available
  oops::Parameter<std::string> ModName{"ModName", "RTTOV", this};

  /// Get the settings for the surface emissivity.  This specifies where
  /// the initial values come from and details about the retrieval.
  oops::Parameter<SurfaceEmissivityParameters> SurfaceEmissivity{
      "surface emissivity", SurfaceEmissivityParameters(), this};

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

  /// Output the Ice Water Path if the profile converges
  oops::Parameter<bool> Store1DVarIWP{"Store1DVarIWP", false, this};

  /// Output the CLW if the profile converges
  oops::Parameter<bool> Store1DVarCLW{"Store1DVarCLW", false, this};

  /// Output the surface to space transmittance if the profile converges
  oops::Parameter<bool> Store1DVarTransmittance{"Store1DVarTransmittance", false, this};

  /// Recalculate the brightness temperature using retrieval values if the
  /// profile has converged
  oops::Parameter<bool> RecalculateBT{"RecalculateBT", false, this};

  /// Flag to read the initial skin temperature from the obs space replacing the
  /// value in the GeoVaLs. The value will then be passed around and updated in the
  /// geovals as the RTTOV interface expects.
  oops::Parameter<bool> SetInitialSkinTFromObsSpace{
      "set the initial skin temperature from the obsspace", false, this};

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

  /// If the iteration number is greater than ConvergeCheckChansAfterIteration then the
  /// slow converging channels, specified by ConvergeCheckChans, have the observation error
  /// inflated to 100000.0
  oops::Parameter<int> ConvergeCheckChansAfterIteration{
      "ConvergeCheckChansAfterIteration", 3, this};

  /// List of channels to inflate the observation error (R) for if the retrieval goes beyond
  /// ConvergeCheckChansAfterIteration iterations.  The inflated variance for these channels is
  /// set to 100000.0 for future iterations effectively removing it from the minimization.
  oops::OptionalParameter<std::vector<int>> ConvergeCheckChans{"ConvergeCheckChans", this};

  /// Check all the retrieved brightness temperatures are within a factor * error of the
  /// observed and bias corrected BTs.  If this value is less than 0.0 this check is
  /// not performed
  oops::Parameter<double> RetrievedErrorFactor{"RetrievedErrorFactor", 4.0, this};

  /// Convergence factor used when the absolute difference in the profile is used
  /// to determine convergence.
  oops::Parameter<double> ConvergenceFactor{"ConvergenceFactor", 0.4, this};

  /// Cost threshold for convergence check when cost function value is used for convergence
  oops::Parameter<double> CostConvergenceFactor{"CostConvergenceFactor", 0.01, this};

  /// The fraction of the Jacobian that is permitted to be below the cloud_top_pressure for the
  /// IR cloudy channel selection.  The Jacobian is integrated from the toa -> surface and a
  /// maximum of 1 % of the integrated Jacobian is allowed to be below the cloud top.
  oops::Parameter<double> IRCloud_Threshold{"IRCloud_Threshold", 0.01, this};

  /// Value to scale the skin temperature error over land.  If less than zero
  /// no scaling is done hence the default value of -1.0.
  oops::Parameter<double> SkinTempErrorLand{"SkinTempErrorLand", -1.0, this};

  /// Max LWP in check the cloudy iteration in kg/m2
  oops::Parameter<double> maxLWPForCloudyCheck{"MaxLWPForCloudyCheck", 2.0, this};

  /// Max IWP in check the cloudy iteration in kg/m2
  oops::Parameter<double> maxIWPForCloudyCheck{"MaxIWPForCloudyCheck", 2.0, this};

  /// -------------------------------
  /// Variables purely for testing
  /// -------------------------------
  /// Starting observation to run through 1d-var, subsetting for testing
  oops::Parameter<int> StartOb{"StartOb", 0, this};

  /// Final observation to run through 1d-var, subsetting for testing
  oops::Parameter<int> FinishOb{"FinishOb", 0, this};
};

}  // namespace ufo

#endif  // UFO_FILTERS_RTTOVONEDVARCHECK_RTTOVONEDVARCHECKPARAMETERS_H_
