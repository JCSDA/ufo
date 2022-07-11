/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSERRORFACTORTOPORAD_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSERRORFACTORTOPORAD_H_

#include <string>
#include <vector>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

///
/// \brief Options applying to observation error inflation as a function of terrain height,
/// channel, and surface-to-space transmittance
///
class ObsErrorFactorTopoRadParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsErrorFactorTopoRadParameters, Parameters)

 public:
  /// List of channels to which the observation error factor applies
  oops::RequiredParameter<std::string> channelList{"channels", this};

  /// Name of the sensor for which the observation error factor applies
  oops::RequiredParameter<std::string> sensor{"sensor", this};

  /// Name of the data group to which the observation error is applied (default: ObsErrorData)
  oops::Parameter<std::string> testObserr{"test_obserr", "ObsErrorData", this};

  /// Name of the data group to which the QC flag is applied  (default is QCflagsData)
  oops::Parameter<std::string> testQCflag{"test_qcflag", "QCflagsData", this};
};

///
/// \brief Error Inflation Factor (EIF) as a function of terrain height, channel,
/// and surface-to-space transmittance
/// H = surface height [m]
/// X = surface-to-space transmittance
/// IASI:
//           EIF = SQRT [ 1 / ( 1 - (1 - (2000/H)^4) * X ] for H > 2000
/// AMSU-A:
///          EIF = SQRT [ 1 / ( 2000 / H ) ] for 2000 < H < 4000 and Channels 1-6,15
///          EIF = SQRT [ 1 / ( 4000 / H ) ] for H > 4000 and Channel 7
/// MHS:
//           EIF = SQRT [ 1 / ( 2000 / H ) ] for H > 2000
///
class ObsErrorFactorTopoRad : public ObsFunctionBase<float> {
 public:
  explicit ObsErrorFactorTopoRad(const eckit::LocalConfiguration &);
  ~ObsErrorFactorTopoRad();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ObsErrorFactorTopoRadParameters options_;
  ufo::Variables invars_;
  std::vector<int> channels_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSERRORFACTORTOPORAD_H_
