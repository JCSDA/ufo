/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_CLWMATCHINDEXMW_H_
#define UFO_FILTERS_OBSFUNCTIONS_CLWMATCHINDEXMW_H_

#include <string>
#include <vector>

#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

class ObsFilterData;

///
/// \brief Options applying to the determination of the cloud match index based on
/// retrieved CLW from observation and background
///
class CLWMatchIndexMWParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(CLWMatchIndexMWParameters, Parameters)

 public:
  /// List of channels available for assimilation
  oops::RequiredParameter<std::string> channelList{"channels", this};

  /// Symmetric cloud amount threshold for each channel
  /// Channel is considered insensitivity to the cloud amount less than the threshold
  oops::RequiredParameter<std::vector<float>> clwretClearSky{"clwret_clearsky", this};

  /// Function to retrieve the cloud liquid water from observation
  oops::RequiredParameter<Variable> clwobsFunction{"clwobs_function", this};

  /// Function to retrieve the cloud liquid water from the simulated observation
  oops::RequiredParameter<Variable> clwbkgFunction{"clwbkg_function", this};
};

///
/// \brief Determination of cloud match index based on retrieved CLW from
/// observation and simulated observation from background
/// 1: both background and observation cont contain clouds
/// 0: either background has cloud or observation has cloud detected
///
class CLWMatchIndexMW : public ObsFunctionBase<float> {
 public:
  explicit CLWMatchIndexMW(const eckit::LocalConfiguration &
                                       = eckit::LocalConfiguration());
  ~CLWMatchIndexMW();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  std::vector<int> channels_;
  CLWMatchIndexMWParameters options_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_CLWMATCHINDEXMW_H_
