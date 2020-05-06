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

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/CLWRetMW.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

///
/// \brief Options applying to the determination of cloud match index based on
/// retrieved CLW from observation and background
///
class CLWMatchIndexMWParameters : public oops::Parameters {
 public:
  /// List of channels available for assimilation
  oops::RequiredParameter<std::string> channelList{"channels", this};

  /// Observation error for each channel under clear-sky condition
  oops::RequiredParameter<std::vector<float>> clwretClearSky{"clwret_clearsky", this};

  /// Function used to retrieve the cloud liquid water from observation (CLWRetMW with ObsValue)
  oops::RequiredParameter<Variable> clwobsFunction{"clwobs_function", this};

  /// Function used to retrieve the cloud liquid water from simulated observation
  /// (CLWRetMW with HofX)
  oops::RequiredParameter<Variable> clwbkgFunction{"clwbkg_function", this};
};

///
/// \brief Determination of cloud match index based on retrieved CLW from
/// observation and simulated observation from background
/// 1: both background and observation cont contain clouds
/// 0: either background has cloud or observation has cloud detected
///
class CLWMatchIndexMW : public ObsFunctionBase {
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
