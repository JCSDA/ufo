/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_CHANNELUSEFLAGCHECKRAD_H_
#define UFO_FILTERS_OBSFUNCTIONS_CHANNELUSEFLAGCHECKRAD_H_

#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

class ObsFilterData;

///
/// \brief Options applying to channel useflag check
///
class ChannelUseflagCheckRadParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ChannelUseflagCheckRadParameters, Parameters)

 public:
  /// List of channels available for assimilation
  oops::RequiredParameter<std::string> channelList{"channels", this};

  /// Useflag (-1: not used; 0: monitoring; 1: used) for each channel in channelList
  oops::RequiredParameter<std::vector<int>> useflagChannel{"use_flag", this};

  /// Configure passive bias correction
  oops::OptionalParameter<bool> passiveBC{"use passive_bc", this};
};

///
/// \brief Channel useflag check: remove channel if useflag is less than one
///
class ChannelUseflagCheckRad : public ObsFunctionBase<float> {
 public:
  explicit ChannelUseflagCheckRad(const eckit::LocalConfiguration &);
  ~ChannelUseflagCheckRad();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ChannelUseflagCheckRadParameters options_;
  ufo::Variables invars_;
  std::vector<int> channels_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_CHANNELUSEFLAGCHECKRAD_H_
