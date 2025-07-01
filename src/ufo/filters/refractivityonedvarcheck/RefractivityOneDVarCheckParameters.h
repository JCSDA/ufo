/*
 * (C) Copyright 2025 UK Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_REFRACTIVITYONEDVARCHECK_REFRACTIVITYONEDVARCHECKPARAMETERS_H_
#define UFO_FILTERS_REFRACTIVITYONEDVARCHECK_REFRACTIVITYONEDVARCHECKPARAMETERS_H_

#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace ufo {

/// Configuration options recognized by the refractivity 1DVar check.
class RefractivityOneDVarCheckParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(RefractivityOneDVarCheckParameters, FilterParametersBase)

 public:
  /// Filename of the static B-matrix to be used in the 1D-Var minimisation
  oops::RequiredParameter<std::string> bmatrix_filename{"bmatrix_filename", this};

  /// Filename of the static R-matrix to be used in the 1D-Var minimisation
  oops::RequiredParameter<std::string> rmatrix_filename{"rmatrix_filename", this};

  /// If true calculate saturation vapour pressure with respect to water and ice
  /// (below zero degrees), else calculate it with respect to water everywhere.
  oops::Parameter<bool> capsupersat{"capsupersat", true, this};

  /// The profile is rejected if the final cost-function is larger than
  /// cost_funct_test times the number of observations.
  oops::Parameter<float> cost_funct_test{"cost_funct_test", 0.5, this};

  /// The minimisation is considered to have converged if the absolute value of
  /// the change in the solution for an iteration, divided by the gradient in
  /// the cost-function is less than Delta_ct2 times the number of observations
  /// in the profile, divided by 200.
  oops::Parameter<float> Delta_ct2{"Delta_ct2", 1, this};

  /// The minimisation is considered to have converged if the absolute change
  /// in the cost-function at this iteration is less than Delta_factor times
  /// either the previous value of the cost-function or the number of
  /// observations (whichever is the smaller). That is, the minimisation has
  /// converged if: ABS(J_new - J_old) < Delta_factor * min(J_old, nObs)
  oops::Parameter<float> Delta_factor{"Delta_factor", 0.01, this};

  /// Threshold for the minimum temperature gradient before a profile is
  /// considered isothermal (units: K per model level). Only applies if
  /// pseudo-levels are being used.
  oops::Parameter<float> min_temp_grad{"min_temp_grad", 1.0e-6, this};

  /// The maximum number of iterations - the profile is rejected if it does not
  /// converge in time.
  oops::Parameter<int> n_iteration_test{"n_iteration_test", 20, this};

  /// If the RMS difference between the observations and the background bending
  /// angle is greater than OB_test then the whole profile is rejected.
  oops::Parameter<float> OB_test{"OB_test", 1, this};

  /// Whether to use pseudo-levels to reduce interpolation errors in the
  /// forward model.
  oops::Parameter<bool> pseudo_ops{"pseudo_ops", true, this};

  /// If true use linear interpolation in ln(pressure), otherwise use linear
  /// interpolation in exner.
  oops::Parameter<bool> vert_interp_ops{"vert_interp_ops", true, this};

  /// If an observation is more than y_test times the observation error away
  /// from the solution bending angle, then the observation (not the whole
  /// profile) is rejected.
  oops::Parameter<float> y_test{"y_test", 5, this};

  /// Minimum threshold for the value of y_test times of observation uncertainty.
  /// Mainly applies to observations at high altitudes where the refractivity is small.
  oops::Parameter<float> minval_ytest{"minval_ytest", 1, this};

  /// Maximum threshold for the value of y_test times of observation uncertainty.
  /// Mainly applies to observations in the lower troposphere where the refractivity is large.
  oops::Parameter<float> maxval_ytest{"maxval_ytest", 30, this};
};

}  // namespace ufo
#endif  // UFO_FILTERS_REFRACTIVITYONEDVARCHECK_REFRACTIVITYONEDVARCHECKPARAMETERS_H_
