/*
 * (C) Crown Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_CLOUDCOSTFUNCTION_H_
#define UFO_FILTERS_OBSFUNCTIONS_CLOUDCOSTFUNCTION_H_

#include <set>
#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

class ObsFilterData;

///
/// \brief Options for calculating Bayesian cost function.
///
class CloudCostFunctionParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(CloudCostFunctionParameters, Parameters)

 public:
  /// Set of channels used in the calculation of the cost function
  oops::RequiredParameter<std::string> chanlist{"cost channels list", this};

  /// Path to location of file describing the R-matrix
  oops::RequiredParameter<std::string> rmatrix_filepath{"RMatrix", this};

  /// Path to location of file describing the B-matrix
  oops::RequiredParameter<std::string> bmatrix_filepath{"BMatrix", this};

  /// List of geovals describing fields required from the B-matrix
  oops::RequiredParameter<std::vector<std::string>>
                                 field_names{"background fields", this};

  /// \brief B-matrix file contains error covariances for ln(qtotal in units g/kg)
  ///
  /// Setting this flag for qtotal requires that the following are all present
  /// in the parameter list "background fields":
  /// - specific_humidity
  /// - mass_content_of_cloud_liquid_water_in_atmosphere_layer
  /// - mass_content_of_cloud_ice_in_atmosphere_layer
  oops::Parameter<bool> qtotal_lnq_gkg{"qtotal", false, this};

  /// Include treatment of rain when splitting total humidity into constituent phases
  oops::Parameter<bool> qtotal_split_rain{"qtotal split rain", false, this};

  /// Include gradient due to ice in brightness temperature total humidity Jacobian
  oops::Parameter<bool> scattering_switch{"scattering radiative transfer", false, this};

  /// Limit specific humidity to minimum value
  oops::Parameter<float> min_q{"minimum specific humidity", 3.0e-6f, this};

  /// Jacobian vertical ordering is reverse of geovals
  oops::Parameter<bool> reverse_Jacobian{"reverse Jacobian order", false, this};

  /// Minimum bound for ObsValue brightness temperature
  oops::Parameter<float> minTb{"minimum ObsValue", 70.0, this};

  /// Maximum bound for ObsValue brightness temperature
  oops::Parameter<float> maxTb{"maximum ObsValue", 340.0, this};

  /// Maximum value for final cost returned by the ObsFunction
  oops::Parameter<float> maxCost{"maximum final cost", 1600.0, this};

  /// \brief Name of the H(x) group used in the cost function calculation.
  ///
  /// H(x) is assumed to be already bias corrected, the default is "HofX"
  ///
  /// Example: use
  ///
  ///          HofX group: MetOfficeBiasCorrHofX
  oops::Parameter<std::string> HofXGroup{"HofX group", "HofX", this};

  /// Vector of channels mapped to surface emissivity error covariances in the B-matrix
  oops::OptionalParameter<std::string> emissMap{"background emissivity channels", this};

  /// Scale B-matrix skin temperature error covariances to user value
  oops::OptionalParameter<float> skinTempError{"skin temperature error", this};
};

///
/// \brief Bayesian cost function for detecting cloud-affected radiances.
///
/// The cloud cost, Jc, is calculated from observation-H(x) departures, y, via
///
/// Jc = (0.5/Nchan) * y.W.y^T
///
/// where Nchan is the number of channels in the calculation and
/// W is the inverse of (H.B.H^T + R):
///
/// H is the Jacobian matrix;
/// B is a background error covariance matrix;
/// R is an observation error covariance matrix.
///
/// The heritage of this code is the Met Office routine Ops_SatRad_CloudCost.
///
/// Implementation here follows Met Office usage, with a static (latitude-varying)
/// B-matrix and a fixed, diagonal R-matrix.
///
/// Reference: S.J. English, J.R. Eyre and J.A. Smith.
/// A cloud‐detection scheme for use with satellite sounding radiances in the
/// context of data assimilation for numerical weather prediction support of
/// nowcasting applications.
/// Quart. J. Royal Meterol. Soc., Vol. 125, pp. 2359-2378 (1999).
/// https://doi.org/10.1002/qj.49712555902

class CloudCostFunction : public ObsFunctionBase<float> {
 public:
  explicit CloudCostFunction(const eckit::LocalConfiguration &
                                       = eckit::LocalConfiguration());
  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;

 private:
  ufo::Variables invars_;
  std::vector<int> channels_;
  std::vector<std::string> fields_;
  std::set<int> emissMap_;
  CloudCostFunctionParameters options_;
};


// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_CLOUDCOSTFUNCTION_H_
