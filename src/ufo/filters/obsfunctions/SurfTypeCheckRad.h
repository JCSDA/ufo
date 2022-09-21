/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_SURFTYPECHECKRAD_H_
#define UFO_FILTERS_OBSFUNCTIONS_SURFTYPECHECKRAD_H_

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
/// \brief Options applying to reject the radiance observations over different surface types
///
class SurfTypeCheckRadParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(SurfTypeCheckRadParameters, Parameters)

 public:
  /// List of channels available for assimilation
  oops::RequiredParameter<std::string> channelList{"channels", this};

  /// Useflag (-1: not used; 0: monitoring; 1: used) for each channel in channelList
  oops::RequiredParameter<std::vector<int>> useflagChannel{"use_flag", this};

  /// convert this number into icld_det,iland_det,isnow_det,imix_det,iice_det,iwater_det...
  oops::RequiredParameter<std::vector<int>> useflagCloudDetect{"use_flag_clddet", this};

  /// Name of the data group to which the observation error is applied (default: ObsErrorData)
  oops::Parameter<std::string> testObserr{"test_obserr", "ObsErrorData", this};

  /// Name of the data group to which the QC flag is applied  (default is QCflagsData)
  oops::Parameter<std::string> testQCflag{"test_qcflag", "QCflagsData", this};
};

///
/// Output of this function:
/// 0 = channel passed by this check
/// 1 = channel rejected by this check
///
class SurfTypeCheckRad : public ObsFunctionBase<float> {
 public:
  explicit SurfTypeCheckRad(const eckit::LocalConfiguration &);
  ~SurfTypeCheckRad();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  SurfTypeCheckRadParameters options_;
  ufo::Variables invars_;
  std::vector<int> channels_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_SURFTYPECHECKRAD_H_
