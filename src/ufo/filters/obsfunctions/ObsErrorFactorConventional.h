/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSERRORFACTORCONVENTIONAL_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSERRORFACTORCONVENTIONAL_H_

#include <memory>
#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

/// \brief Options controlling ObsErrorFactorConventional ObsFunction
class ObsErrorFactorConventionalParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsErrorFactorConventionalParameters, Parameters)

 public:
  /// Inflate variables
  oops::RequiredParameter<std::vector<std::string>> inflatevars{"inflate variables", this};
  /// Name of the data group to which the QC flag is applied  (default is QCflagsData)
  oops::Parameter<std::string> testQCflag{"test QCflag", "QCflagsData", this};
  /// Data with QC <= threshold values will be used in this obsFunction.
  /// Only used when "test QCflag" is "PreQC" (include pre-defined QC flags)
  oops::OptionalParameter<int> qcthreshold{"test QCthreshold", this};
  /// Threshold for horizontal distance or chord distance. While the horizontal(chord)
  /// distance between the two vertical points is bigger than the threshold,
  /// these observations will not be considered for the obs inflation
  oops::OptionalParameter<float> distthreshold{"distance threshold", this};
  /// Name of air pressure
  oops::Parameter<std::string> pressureFullName{"pressure",
                                                "Full name of air pressure",
                                                "ObsValue/pressure",
                                                this};
};

// -----------------------------------------------------------------------------

/// \brief Inflate the observation error for conventional as done by GSI-Observer.
///
/// This routine was designed to mimic the GSI observer code (i.e., subroutine errormod in
/// qcmod.f90) to inflate the observation error for conventional/satwinds using the
/// qc flags generated from a filter or from the input files with a group name
/// (test QCflag) defined in the yaml. The inflation factor is determined by the obseravtion
/// vertical spacing (in pressure) relative to the corresponding model pressure interval.
/// This error inflation obsFuction is used in GSI for temperature, moiture, and winds from
/// conventional obs as well as some of satellite retrievels, e.g., radiosonde/other conventional
/// temperature, moisture, and wind, SCAT winds, VAD winds, and potentially aircraft ascent
/// and descent profiles.
///
/// Notes:
/// (1) If using this obs function in a filter, please
/// make sure the "filter variables" and "inflate variables" prescribed
/// with the same variable name.
/// (2) This obs funciton requires each of the obs profiles are sorted by pressure
/// in descending order
///
/// ### example configurations for testing this obs function: ###
///
///  obs function:
///    name: ObsFunction/ObsErrorFactorConventional
///    variables: [windEastward]   # Variable name for output
///    tolerance: 1.e-6
///    options:
///      inflate variables: [windEastward]  # Ok to be multiple dimensions for running
///                                         # this obsFunction only (not within a filter)
///      test QCflag: PreQC  # Optional. If not defined, use QCflags from prior filters
///      test QCthreshold: 2 # Optonal, only when PreQC is used
///                          # Data with PreQC flags <= 2 will be used in this ObsFunction
///                          # Default is 3 for PreQC
///                          # In GSI(PreQC): if noiqc (no oiqc)=true, QCthreshold=7;
///                                           if noiqc=false, QCthreshold=3
///
/// ### example configurations for using this obs function in a filter: ###
///   - filter: BlackList
///     filter variables:
///     - name: virtualTemperature # Have to be consistent with "inflate
///                                 # variables". Therefore, only one variable allowed
///                                 # while running with this obsFunc
///     action:
///       name: inflate error
///       inflation variable:
///         name: ObsFunction/ObsErrorFactorConventional
///         options:
///           inflate variables: [virtualTemperature]  # Have to be consistent with "filter
///                                                     # variables". Therefore, only one
///                                                     # variable allowed
///
/// ### example configurations for using obsgrouping: ###
///
///      obsgrouping:
///        group variables: ["stationIdentification", "dateTime"] # Choose parameteres to identify
///                                                    # each of the obs profiles
///        sort variable: "pressure"
///        sort order: "descending"
///
class ObsErrorFactorConventional : public ObsFunctionBase<float> {
 public:
  static const std::string classname() {return "ObsErrorFactorConventional";}

  explicit ObsErrorFactorConventional(const eckit::Configuration &config);
  ~ObsErrorFactorConventional();

  void compute(const ObsFilterData &, ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  std::unique_ptr<ObsErrorFactorConventionalParameters> options_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSERRORFACTORCONVENTIONAL_H_
