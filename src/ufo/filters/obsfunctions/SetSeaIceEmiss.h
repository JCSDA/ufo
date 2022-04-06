/*
 * (C) Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_SETSEAICEEMISS_H_
#define UFO_FILTERS_OBSFUNCTIONS_SETSEAICEEMISS_H_

#include <string>
#include <vector>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/utils/SurfaceReportConstants.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

///
/// \brief Options applying to the setting of sea-ice emissivity
///
class SetSeaIceEmissParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SetSeaIceEmissParameters, Parameters)

 public:
  /// List of channels
  oops::RequiredParameter<std::string> channelList{"channels", this};

  /// Polarization indices (0[0.5(V+H)],1[V],2[H]) of each channels
  oops::RequiredParameter<std::vector<size_t>> polList{"polarization index", this};

  /// Central frequency for channel in GHz
  oops::RequiredParameter<std::vector<float>> freqList{"channel frequency", this};

  /// Nominal orbit height
  oops::RequiredParameter<float> orbitHeight{"orbit height", this};
};

///
/// \brief Set emissivity according to FASTEM for sea-ice scenes based on model and observation data
///
class SetSeaIceEmiss : public ObsFunctionBase<float> {
 public:
  explicit SetSeaIceEmiss(const eckit::LocalConfiguration &);
  ~SetSeaIceEmiss();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;

  const ufo::Variables & requiredVariables() const;

 private:
  ufo::Variables invars_;
  std::vector<int> channels_;
  SetSeaIceEmissParameters options_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_SETSEAICEEMISS_H_
