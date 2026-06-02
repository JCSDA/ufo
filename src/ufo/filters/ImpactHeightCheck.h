/*
 * (C) British Crown Copyright 2021 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_IMPACTHEIGHTCHECK_H_
#define UFO_FILTERS_IMPACTHEIGHTCHECK_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/NumericConstraints.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "oops/util/Printable.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

/// Parameters controlling the operation of the ImpactHeightCheck filter.
class ImpactHeightCheckParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(ImpactHeightCheckParameters, FilterParametersBase)

 public:
  /// The threshold used to define a sharp gradient in refractivity.
  /// Units: N-units / m
  /// Inversions, and possible ducting, are identified by looking for sharp gradients
  /// in the refractivity.  If the refractivity gradient is less than this, then any
  /// data below this point are rejected.
  oops::Parameter<float> gradientThreshold{"gradient threshold", -0.08f, this};

  /// The height (in m) of a buffer-zone for rejecting data above sharp gradients.
  /// If a sharp gradient in refractivity is identified, then all data above
  /// this point (plus sharpGradientOffset) is rejected.  This parameter makes
  /// sure that we don't use any data which is too close to the sharp gradient.
  oops::Parameter<float> sharpGradientOffset{"sharp gradient offset", 500, this};

  /// Reject data within this height (in m) of the surface.
  oops::Parameter<float> surfaceOffset{"surface offset", 600, this};

  /// The maximum impact height (in m) at which to accept observations
  oops::Parameter<float> maximumHeight{"maximum height", 80000, this};

  /// Whether to print extra verbose output from this routine.
  oops::Parameter<bool> verboseOutput{"verbose output", false, this};
};

/// ImpactHeightCheck: Calculate the impact height for model profiles, and reject
/// any observations which are outside the range of model impact heights.  Check
/// for any sharp refractivity gradients, and reject any observations below them.
/// Refractivity and model heights must have been saved into ObsDiagnostics by
/// the observation operator.

class ImpactHeightCheck : public FilterBase,
                           private util::ObjectCounter<ImpactHeightCheck> {
 public:
  /// The type of parameters accepted by the constructor of this filter.
  /// This typedef is used by the FilterFactory.
  typedef ImpactHeightCheckParameters Parameters_;

  static const std::string classname() {return "ufo::ImpactHeightCheck";}

  ImpactHeightCheck(ioda::ObsSpace &, const Parameters_ &,
                     std::shared_ptr<ioda::ObsDataVector<int> >,
                     std::shared_ptr<ioda::ObsDataVector<float> >);
  ~ImpactHeightCheck();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::domain;}
  Parameters_ parameters_;
  std::vector<float> calcVerticalGradient(const std::vector<float> &,
                                          const std::vector<float> &) const;
  float calcImpactHeight(float, float, float) const;
};

}  // namespace ufo

#endif  // UFO_FILTERS_IMPACTHEIGHTCHECK_H_
