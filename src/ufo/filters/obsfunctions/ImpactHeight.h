/* -----------------------------------------------------------------------------
 * (C) British Crown Copyright 2020 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * -----------------------------------------------------------------------------
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_IMPACTHEIGHT_H_
#define UFO_FILTERS_OBSFUNCTIONS_IMPACTHEIGHT_H_

#include <string>
#include <vector>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {
///
/// \brief Options: list of channels to apply to
///
class ImpactHeightParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ImpactHeightParameters, Parameters)

 public:
  /// List of channels available for assimilation
  oops::Parameter<std::string> channelList{"channels", "", this};
};

///
/// \brief Function calculates the GNSS-RO impact height as the difference
/// between the impact parameter and earth's radius of curvature.
///
class ImpactHeight : public ObsFunctionBase<float> {
 public:
  explicit ImpactHeight(const eckit::LocalConfiguration &);
  ~ImpactHeight();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  std::vector<int> channels_;
  ImpactHeightParameters options_;
};
// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_IMPACTHEIGHT_H_
