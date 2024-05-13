/*
 * (C) Copyright 2017-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 * This filter implements the super refraction check used in Meteo France
 * Created following the advice of Dominique Raspaud (Meteo France) and Neill Bowler (Met Office)
 * This method checks the observed refractivity profiles and accepts an observation only if
 * all the following criteria are met: dN/dz>-50 km^(-1), dN/dz<-0.001 km^(-1), dN/dz≠0 and |(d^2 N)/dz^2 |<100 km^(-2)
 */

#ifndef UFO_FILTERS_OBSREFRACTIVITYGRADIENTCHECK_H_
#define UFO_FILTERS_OBSREFRACTIVITYGRADIENTCHECK_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"
namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

/// GNSSRO super refraction check based on observed refractivity gradient
class ObsRefractivityGradientCheckParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsRefractivityGradientCheckParameters, FilterParametersBase)

 public:
  oops::Parameter<float> gradientMin{"gradient min",
                      "The minimum allowed refractivity gradient. "
                      "Gradients below this value are marked as having super-refraction.",
                       -0.05, this};
  oops::Parameter<float> gradientMax{"gradient max",
                      "The maximum allowed refractivity gradient. "
                      "Above this value the refractivity gradient is not changing quickly enough, "
                      "which can be associated with super-refraction.",
                      -1.0e-6, this};
  oops::Parameter<float> secondDerivative{"second derivative",
                      "The maximum value for the second derivative check. "
                       "If the absolute value of the second derivative is above this, "
                      "then the observation will be rejected. ",
                       0.0001, this};
  oops::Parameter<float> maxCheckHeight{"max check height",
                       "The height below what the check is applied (unit: meter). ",
                        18000, this};
};

class ObsRefractivityGradientCheck : public FilterBase,
      private util::ObjectCounter<ObsRefractivityGradientCheck> {
 public:
  typedef ObsRefractivityGradientCheckParameters Parameters_;
  static const std::string classname() {return "ufo::ObsRefractivityGradientCheck";}

  ObsRefractivityGradientCheck(ioda::ObsSpace &, const Parameters_ &,
                  std::shared_ptr<ioda::ObsDataVector<int> >,
                  std::shared_ptr<ioda::ObsDataVector<float> >);
  ~ObsRefractivityGradientCheck();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::superrefraction;}
  std::vector<float> calcVerticalGradient(const std::vector<float> &,
                                          const std::vector<float> &) const;
  Parameters_ parameters_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_OBSREFRACTIVITYGRADIENTCHECK_H_
