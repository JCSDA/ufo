/*
 * (C) Copyright 2023 NOAA/NWS/NCEP/EMC
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_AIRCRAFTBIASCORRECTIONTERM_H_
#define UFO_FILTERS_OBSFUNCTIONS_AIRCRAFTBIASCORRECTIONTERM_H_

#include <string>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

///
/// \brief This function requires the following required parameters:
///        coeff_grpvarname (std::string): <Group>/<Variable> of BC coefficient
///        predi_grpvarname (std::string): <Group>/<Variable> of BC predictor
///        predi_order (std::float): order to which BC predictor is raised
///
class AircraftBiasCorrectionTermParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(AircraftBiasCorrectionTermParameters, Parameters)

 public:
  // <Group>/<Variable> of BC coefficient
  oops::RequiredParameter<std::string> coeff_grpvarname{"coeff_grpvarname", this};
  // <Group>/<Variable> of BC predictor
  oops::RequiredParameter<std::string> predi_grpvarname{"predi_grpvarname", this};
  // Order to which BC predictor is raised
  oops::RequiredParameter<float> predi_order{"predi_order", this};
};

// -----------------------------------------------------------------------------

/// \brief Compute the bias correction term for aircraft bias correction, which is
///        equal to <coefficient> * std::pow(<predictor>,<order>). This can be
///        combined with the aircraft airTemperature and other predictors through
///        use of a LinearCombination filter to bias correct the airTemperature
///        observations.
///
/// ~~~
///
/// ### Sample YAML configuration: Generate BC term for aircraft ascent term (order:1)
///     - filter: Variable Assignment
///       filter variables:
///       - name: airTemperature
///       assignments:
///       - name: BiasCorrectionTerm/ascentPredictor
///         type: float
///       function:
///         name: ObsFunction/AircraftBiasCorrectionTerm
///         options:
///           coeff_grpvarname: BiasCorrectionCoefficient/ascentPredictor
///           predi_grpvarname: MetaData/windUpward
///           predi_order: 1.
///
class AircraftBiasCorrectionTerm : public ObsFunctionBase<float> {
 public:
  static const std::string classname() {return "AircraftBiasCorrectionTerm";}

  explicit AircraftBiasCorrectionTerm(const eckit::LocalConfiguration &
                                 = eckit::LocalConfiguration());
  ~AircraftBiasCorrectionTerm();

  void compute(const ObsFilterData &, ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  AircraftBiasCorrectionTermParameters options_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_AIRCRAFTBIASCORRECTIONTERM_H_
