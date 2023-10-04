/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_NEARSSTRETCHECKIR_H_
#define UFO_FILTERS_OBSFUNCTIONS_NEARSSTRETCHECKIR_H_

#include <string>
#include <vector>

#include "boost/optional/optional_io.hpp"

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

///
/// \brief Options controlling the quality control steps using retrieved near-sea-surface
/// temperature (NSST) for Infrared sensors
///
class NearSSTRetCheckIRParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(NearSSTRetCheckIRParameters, Parameters)

 public:
  /// List of channels available for assimilation
  oops::RequiredParameter<std::string> channelList{"channels", this};

  /// Useflag (-1: not used; 0: monitoring; 1: used) for each channel in channelList
  oops::RequiredParameter<std::vector<int>> useflagChannel{"use_flag", this};

  /// Name of the data group to which the observation error is applied (default: ObsErrorData)
  oops::Parameter<std::string> testObserr{"test_obserr", "ObsErrorData", this};

  /// Name of the HofX group used to replace the default group (default is HofX)
  oops::Parameter<std::string> testHofX{"test_hofx", "HofX", this};

  /// Name of the data group to which the QC flag is applied  (default is QCflagsData)
  oops::Parameter<std::string> testQCflag{"test_qcflag", "QCflagsData", this};

  /// Parameter for original observation error
  oops::OptionalParameter<std::vector<float>> obserrOriginal{"error parameter vector", this};
};

///
/// \brief QC using retrieved near-sea-surface temperature (NSST) from Infrared radiances
/// 2-step QC procedures:
/// (1) Perform NSST retrieval from radiances at obs location, and obtained
///     increment of NSST from its first guess value
/// (2) For surface sensitive channels, remove them from assimilation if the
///     retrieved NSST increment from step (1) is larger than a pre-defined
///     threshold
/// Output:
/// 0 = channel is retained for assimilation
/// 1 = channel is not retained for assimilation
///
class NearSSTRetCheckIR : public ObsFunctionBase<float> {
 public:
  explicit NearSSTRetCheckIR(const eckit::LocalConfiguration &);
  ~NearSSTRetCheckIR();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  NearSSTRetCheckIRParameters options_;
  ufo::Variables invars_;
  std::vector<int> channels_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_NEARSSTRETCHECKIR_H_
