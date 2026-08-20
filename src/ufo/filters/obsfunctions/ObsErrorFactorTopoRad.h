/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSERRORFACTORTOPORAD_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSERRORFACTORTOPORAD_H_

#include <map>
#include <string>
#include <vector>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

///
/// \brief A group of channels that share the same terrain-height reference (threshold) used
/// by ObsErrorFactorTopoRad. Channels not covered by any group are left uncorrected.
///
class ObsErrorFactorTopoRadHeightGroupParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsErrorFactorTopoRadHeightGroupParameters, Parameters)

 public:
  /// List of channels (e.g. "1-5, 9, 12, 15, 20-24") to which this reference height applies
  oops::RequiredParameter<std::string> channelList{"channels", this};

  /// Terrain-height reference in meters for these channels (default: 2000)
  oops::Parameter<float> height{"height", 2000.0, this};
};

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

  /// Optional list of channel/reference-height groups, e.g.:
  ///   height groups:
  ///     - channels: 1-5, 6, 9
  ///       height: 2000
  ///     - channels: 8
  ///       height: 4000
  /// Channels not listed in any group receive no correction (error inflation factor of 1).
  /// If this list is left empty (the default), the correction is applied to every channel in
  /// `channels` using a reference height of 2000 m.
  oops::Parameter<std::vector<ObsErrorFactorTopoRadHeightGroupParameters>>
       heightGroups{"height groups", {}, this};
};

///
/// \brief Error Inflation Factor (EIF) as a function of terrain height, channel,
/// and surface-to-space transmittance
/// H = surface height [m]
/// X = surface-to-space transmittance
/// R = terrain-height reference (threshold) for the channel [m] (see `height groups`, default
///     2000 for every channel)
/// Infrared:
//           EIF = SQRT [ 1 / ( 1 - (1 - (R/H)^4) * X ] for H > R
/// Microwave:
///          EIF = SQRT [ 1 / ( R / H ) ] for H > R
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
  std::string sensor_;
  /// Maps each channel that should be corrected to its terrain-height reference (threshold).
  /// Channels absent from this map are left uncorrected.
  std::map<int, float> channelHeights_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSERRORFACTORTOPORAD_H_
