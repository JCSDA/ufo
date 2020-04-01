/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONCLWRETMEAN_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONCLWRETMEAN_H_

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
/// \brief Calculate symmetric (mean) cloud amount from the cloud amount retrieved
/// from the observed and simulated measurements
///
class ObsFunctionCLWRetMeanParameters : public oops::Parameters {
 public:
  ///
  /// Required Parameters:
  ///
  /// Names of the data group used to retrieve the cloud liquid water
  /// Example: get retrieved CLW from observation and simulated observation respectively
  ///          clwret_types: [ObsValue, HofX]
  /// Example: get retrieved CLW from observation or simulated observation only
  ///          clwret_types: [ObsValue]
  ///          clwret_types: [HofX]
  oops::RequiredParameter<std::vector<std::string>> varGrp{"clwret_types", this};
  ///
};

class ObsFunctionCLWRetMean : public ObsFunctionBase {
 public:
  explicit ObsFunctionCLWRetMean(const eckit::LocalConfiguration &
                                       = eckit::LocalConfiguration());
  ~ObsFunctionCLWRetMean();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  ObsFunctionCLWRetMeanParameters options_;
  eckit::LocalConfiguration conf_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONCLWRETMEAN_H_
