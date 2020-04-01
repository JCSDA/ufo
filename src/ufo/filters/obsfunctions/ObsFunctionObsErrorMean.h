/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONOBSERRORMEAN_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONOBSERRORMEAN_H_

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
/// \brief Calculate symmetric (mean) observation error from the average of the observed and
/// simulated cloud amount.
///
/// References:
/// (1) Geer, A. J. and P. Bauer (2011). Observation errors in all-sky data
///     assimilation. Quart. J. Roy. Meteorol. Soc. 137, 2024–2037.
/// (2) Geer, A. J., P. Bauer, and S. J. English (2012). Assimilating AMSU-A
///     temperature sounding channels in the presence of cloud and precipitation.
///     Published simultaneously as ECMWF Technical Memoranda 670 and ECMWF/EUMETSAT
///     fellowship reports 24.
///
class ObsFunctionObsErrorMeanParameters : public oops::Parameters {
 public:
  ///
  /// Required Parameters:
  ///
  /// List of channels to which the observation error estimation applies
  oops::RequiredParameter<std::string> chList{"channels", this};
  ///
  /// channel number corresponding to 23.8GHz to which the retrieval
  /// of cloud liquid water applies
  /// Example: AMSU-A channel numbers used in cloud liquid water retrieval
  ///          clwret_channels: 1
  ///
  oops::RequiredParameter<int> ch238{"clwret_ch238", this};
  ///
  /// channel number corresponding to 31.4GHz to which the retrieval
  /// of cloud liquid water applies
  /// Example: AMSU-A channel numbers used in cloud liquid water retrieval
  ///          clwret_channels: 2
  oops::RequiredParameter<int> ch314{"clwret_ch314", this};
  ///
  /// names of the data group used for the retrieval of cloud liquid water
  /// vargrp: [ObsValue, HofX]
  oops::RequiredParameter<std::vector<std::string>> varGrp{"clwret_types", this};
  ///
  /// Symmetric cloud liquid water threshold corresponding to minimum observation error
  oops::RequiredParameter<std::vector<float>> clwClr{"clw_clr_threshold", this};
  ///
  /// Symmetric cloud liquid water threshold correspondding to maximum observation error
  oops::RequiredParameter<std::vector<float>> clwCld{"clw_cld_threshold", this};
  ///
  /// Minimum value for observation error
  oops::RequiredParameter<std::vector<float>> obserrMin{"obserr_min", this};
  ///
  /// Maximum value for observation error
  oops::RequiredParameter<std::vector<float>> obserrMax{"obserr_max", this};
  ///
  /// Optional Parameters :
  ///
  /// Name of the data group to which the bias correction is applied (default: no bias applied)
  /// Example: add bias corretion to simulated observation
  ///          bias_application: HofX
  /// Example: add bias corretion to observation
  ///          bias_application: ObsValue
  oops::Parameter<std::string> addBias{"bias_application", std::string(), this};
  ///
  /// Name of the bias correction group used to replace the default group (default group is ObsBias)
  /// Example: use observation bias correction values from GSI
  ///          test_groups: GsiObsBias
  oops::Parameter<std::string> testGrp{"test_group", std::string(), this};
  ///
};

class ObsFunctionObsErrorMean : public ObsFunctionBase {
 public:
  explicit ObsFunctionObsErrorMean(const eckit::LocalConfiguration &);
  ~ObsFunctionObsErrorMean();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  eckit::LocalConfiguration conf_;
  ObsFunctionObsErrorMeanParameters options_;
  ufo::Variables invars_;
  std::vector<int> channels_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONOBSERRORMEAN_H_
