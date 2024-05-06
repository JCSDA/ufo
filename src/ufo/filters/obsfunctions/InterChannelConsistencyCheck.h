/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_INTERCHANNELCONSISTENCYCHECK_H_
#define UFO_FILTERS_OBSFUNCTIONS_INTERCHANNELCONSISTENCYCHECK_H_

#include <string>
#include <vector>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

///
/// \brief Options applying to inter-channel consistency check
///
class InterChannelConsistencyCheckParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(InterChannelConsistencyCheckParameters, Parameters)

 public:
  /// List of channels available for assimilation
  oops::RequiredParameter<std::string> channelList{"channels", this};

  /// Name of the sensor for which the observation error factor applies
  oops::RequiredParameter<std::string> sensor{"sensor", this};

  /// Useflag (-1: not used; 0: monitoring; 1: used) for each channel in channelList
  oops::RequiredParameter<std::vector<int>> useflagChannel{"use_flag", this};

  /// Configure passive bias correction
  oops::Parameter<bool> passiveBC{"use passive_bc", false, this};

  /// Name of the data group to which the observation error is applied (default: ObsErrorData)
  oops::Parameter<std::string> testObserr{"test_obserr", "ObsErrorData", this};

  /// Name of the data group to which the QC flag is applied  (default is QCflagsData)
  oops::Parameter<std::string> testQCflag{"test_qcflag", "QCflagsData", this};
};

///
/// \brief Inter-channel consistency check
///
class InterChannelConsistencyCheck : public ObsFunctionBase<float> {
 public:
  explicit InterChannelConsistencyCheck(const eckit::LocalConfiguration &);
  ~InterChannelConsistencyCheck();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  std::vector<int> channels_;
  InterChannelConsistencyCheckParameters options_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_INTERCHANNELCONSISTENCYCHECK_H_
