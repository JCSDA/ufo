/*
 * (C) Copyright 2017-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 * This filter implements the super refraction check used in NCEP (GSI V16.3)
 * This method checks if the refractivity grdients between model levels are
 * greater than a series thresholds
 */

#ifndef UFO_FILTERS_SUPERREFRACTIONCHECKNBAM_H_
#define UFO_FILTERS_SUPERREFRACTIONCHECKNBAM_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/NumericConstraints.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"
namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

/// GNSSRO super refraction check implemented in NCEP/EMC (GSI V16.3)
class SuperRefractionCheckNBAMParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(SuperRefractionCheckNBAMParameters, FilterParametersBase)

 public:
  oops::Parameter<int> step{"step",
                      "steps to check the super refraction"
                      "(1-single step, 2-dual, or 3-dual alternative)",
                       2, this};
  oops::Parameter<float> checkHeight{"check height",
                        "below which height (unit: meter) the super refraction check is applied",
                         5000, this};
  oops::Parameter<int> closeLayers{"close layers",
                      "number of model layers close-to-SR which will also be rejected",
                       5, this};
  oops::Parameter<float> thresholdS1{"threshold step1",
                        "Threshold for super-refraction which will be used in step 1"
                        " (fraction of critical value)",
                         0.75, this};
  oops::Parameter<float> thresholdS2{"threshold step2",
                        "Threshold for super-refraction which will be used in step 2"
                        " (fraction of critical value)",
                         0.5, this};
  oops::Parameter<float>maxObs{"max obs value",
                       "maximum obs value to compare in step2 check",
                         0.03, this};
};

class SuperRefractionCheckNBAM : public FilterBase,
      private util::ObjectCounter<SuperRefractionCheckNBAM> {
 public:
  typedef SuperRefractionCheckNBAMParameters Parameters_;
  static std::string classname() {return std::string("ufo::SuperRefractionCheckNBAM");}

  SuperRefractionCheckNBAM(ioda::ObsSpace &, const Parameters_ &,
                  std::shared_ptr<ioda::ObsDataVector<int> >,
                  std::shared_ptr<ioda::ObsDataVector<float> >);
  ~SuperRefractionCheckNBAM();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::superrefraction;}
  std::vector<float> calcImpactParameterModel(const std::vector<float> &,
                                              const std::vector<float> &,
                                              float, float, float) const;
  std::vector<float> calcRefractivityGradient(const std::vector<float> &,
                                              const std::vector<float> &) const;
  Parameters_ parameters_;
};

}  //  namespace ufo

#endif  // UFO_FILTERS_SUPERREFRACTIONCHECKNBAM_H_
