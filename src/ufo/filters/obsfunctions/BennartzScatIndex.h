/*
 * (C) Copyright 2019 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_BENNARTZSCATINDEX_H_
#define UFO_FILTERS_OBSFUNCTIONS_BENNARTZSCATINDEX_H_

#include <string>
#include <vector>

#include "oops/util/parameters/NumericConstraints.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

class ObsFilterData;

///
/// \brief Options for calculating scattering index from 89 GHz and 150 GHz channels.
///
class BennartzScatIndexParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(BennartzScatIndexParameters, Parameters)

 public:
  /// \brief Channel number corresponding to 89 GHz
  /// (or nearby frequency depending on sensor channel specification)
  ///
  /// Example: MHS channel 1 at 89 GHz
  ///
  ///          channel_89ghz: 1
  oops::RequiredParameter<int> ch89{"channel_89ghz", this, {oops::minConstraint(1)}};

  /// \brief Channel number corresponding to 150 GHz
  /// (or nearby frequency depending on sensor channel specification)
  ///
  /// Example: MHS channel 2 at 157 GHz
  ///
  ///          channel_150ghz: 2
  oops::RequiredParameter<int> ch150{"channel_150ghz", this, {oops::minConstraint(1)}};

  /// \brief First coefficient used to compute scattering index
  ///
  /// Offset term is calculated as coeff1 + coeff2*sensor_zenith_angle
  oops::RequiredParameter<float> coeff1{"bennartz_coeff_1", this};

  /// \brief Second coefficient used to compute scattering index
  ///
  /// Offset term is calculated as coeff1 + coeff2*sensor_zenith_angle
  oops::RequiredParameter<float> coeff2{"bennartz_coeff_2", this};

  /// \brief Name of the bias correction group used to apply correction to ObsValue.
  ///
  /// Default (missing optional parameter) applies no bias
  ///
  /// Example: use ObsBias correction values
  ///
  ///          apply_bias: ObsBias
  oops::Parameter<std::string> applyBias{"apply_bias", "", this};
};

///
/// \brief Calculate scattering index from 89 GHz and 150 GHz channels.
///
/// Reference: R. Bennartz (2002)
/// Precipitation analysis using the Advanced Microwave Sounding Unit in
/// support of nowcasting applications
/// Meteorol. Appl. 9, 177-189 (2002) DOI:10.1017/S1350482702002037

class BennartzScatIndex : public ObsFunctionBase<float> {
 public:
  explicit BennartzScatIndex(const eckit::LocalConfiguration &
                                       = eckit::LocalConfiguration());
  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  BennartzScatIndexParameters options_;
  std::vector<int> channels_ = {0, 0};
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_BENNARTZSCATINDEX_H_
