/* -----------------------------------------------------------------------------
 * (C) British Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_SATWINDINVERSIONCORRECTION_H_
#define UFO_FILTERS_SATWINDINVERSIONCORRECTION_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

// Define cloud motion methods
enum CloudMotionMethod {
  infrared    = 1,  /// Motion observed in the infrared channel
  visible     = 2,  /// Motion observed in the visible channel
  vapourcloud = 3,  /// Motion observed in the water vapour channel in cloud
  combination = 4,  /// Motion observed in a combination of spectral channels
  vapourclear = 5,  /// Motion observed in the water vapour channel in clear air
  ozone       = 6,  /// Motion observed in the ozone channel
  vapour      = 7   /// Motion observed in the water vapour channel (cloud or clear)
};

/// \brief Parameters controlling the operation of the SatwindInversionCorrection filter.
class SatwindInversionCorrectionParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(SatwindInversionCorrectionParameters, FilterParametersBase)

 public:
  /// Name of the observation pressure variable to correct
  oops::RequiredParameter<Variable> obs_pressure{"observation pressure", this};
  /// Relative humidity (%) threshold value
  oops::RequiredParameter<float> rh_threshold{"RH threshold", this};
  /// Minimum AMV pressure (Pa) to consider for correction - set default
  oops::Parameter<float> min_pressure{"minimum pressure", 70000.0, this};
  /// Maximum model pressure (Pa) to consider - set default
  oops::Parameter<float> max_pressure{"maximum pressure", 105000.0, this};
  /// Temperature difference (K) between inversion base and top - set default
  oops::Parameter<float> inversion_temperature{"inversion temperature", 2.0, this};
};

// -----------------------------------------------------------------------------

/// \brief A filter that modifies the assigned pressure of AMV observations if a
/// temperature inversion is detected in the model profile and defined criteria are met.
///
/// See SatwindInversionCorrectionParameters for the documentation of the parameters controlling
/// this filter.
class SatwindInversionCorrection : public FilterBase,
                       private util::ObjectCounter<SatwindInversionCorrection> {
 public:
  /// The type of parameters accepted by the constructor of this filter.
  /// This typedef is used by the FilterFactory.
  typedef SatwindInversionCorrectionParameters Parameters_;

  static const std::string classname() {return "ufo::SatwindInversionCorrection";}

  SatwindInversionCorrection(ioda::ObsSpace &, const Parameters_ &,
                  std::shared_ptr<ioda::ObsDataVector<int> >,
                  std::shared_ptr<ioda::ObsDataVector<float> >);
  ~SatwindInversionCorrection();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::pass;}

  Parameters_ parameters_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_SATWINDINVERSIONCORRECTION_H_
