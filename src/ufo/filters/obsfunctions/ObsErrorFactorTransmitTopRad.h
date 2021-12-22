/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSERRORFACTORTRANSMITTOPRAD_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSERRORFACTORTRANSMITTOPRAD_H_

#include <string>
#include <vector>

#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

class ObsFilterData;

///
/// \brief Options applying to observation error inflation as a function to
/// model top-to-space transmittance
///
class ObsErrorFactorTransmitTopRadParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsErrorFactorTransmitTopRadParameters, Parameters)

 public:
  /// List of channels to which the observation error factor applies
  oops::RequiredParameter<std::string> channelList{"channels", this};
};

///
/// \brief Error Inflation Factor (EIF) for satellite radiance as a function of model
/// top-to-space transmittance:
/// tao = model top-to-space transmittance
/// EIF = SQRT ( 1.0 / tao )
///
class ObsErrorFactorTransmitTopRad : public ObsFunctionBase<float> {
 public:
  explicit ObsErrorFactorTransmitTopRad(const eckit::LocalConfiguration &);
  ~ObsErrorFactorTransmitTopRad();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ObsErrorFactorTransmitTopRadParameters options_;
  ufo::Variables invars_;
  std::vector<int> channels_;
  const eckit::LocalConfiguration conf_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSERRORFACTORTRANSMITTOPRAD_H_
