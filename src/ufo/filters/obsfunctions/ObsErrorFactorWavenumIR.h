/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSERRORFACTORWAVENUMIR_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSERRORFACTORWAVENUMIR_H_

#include <string>
#include <vector>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

///
/// \brief Error Inflation Factor (EIF) for satellite radiance with wavenumber in the
/// range of (2000, 2400] during daytime (sun zenith angle < 89) and containing
/// water fraction in the footprint:
/// x = wavenumber [1/cm]
/// y = surface-to-space transmittance
/// z = solar zenith angle [radian]
/// EIF = SQRT[ 1 / ( 1 - (x - 2000)) * y * MAX(0, COS(z)) / 4000 ]
///
class ObsErrorFactorWavenumIRParameters : public oops::Parameters {
 public:
  ///
  /// Required Parameters:
  ///
  /// List of channels to which the observation error factor applies
  oops::RequiredParameter<std::string> chlist{"channels", this};
  ///
};

class ObsErrorFactorWavenumIR : public ObsFunctionBase {
 public:
  explicit ObsErrorFactorWavenumIR(const eckit::LocalConfiguration &);
  ~ObsErrorFactorWavenumIR();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ObsErrorFactorWavenumIRParameters options_;
  ufo::Variables invars_;
  std::vector<int> channels_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSERRORFACTORWAVENUMIR_H_
