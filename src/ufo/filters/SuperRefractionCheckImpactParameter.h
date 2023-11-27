/*
 * (C) Copyright 2017-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 * This filter implements the super refraction check used in ECMWF
 * and Met Office.
 * This method checks if the vertical difference of the modelled impact parameter
 * between a given layer and the one below is larger than the threshold (default is 10m).
 * observations below the identified model height are marked as having super-refraction.
 */

#ifndef UFO_FILTERS_SUPERREFRACTIONCHECKIMPACTPARAMETER_H_
#define UFO_FILTERS_SUPERREFRACTIONCHECKIMPACTPARAMETER_H_

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

/// GNSSRO super refraction check implemented in UKMO and ECMWF
class SuperRefractionCheckImpactParameterParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(SuperRefractionCheckImpactParameterParameters, FilterParametersBase)

 public:
  oops::Parameter<float> threshold{"threshold",
                      "The maximum allowed impact parameter difference in model background. "
                      "Observations above are marked as having super-refraction.",
                       10, this};
  oops::Parameter<bool> profileCheck{"profile check",
                      "if only use the geovals of the bottom observation of the profile",
                       false, this};
};

class SuperRefractionCheckImpactParameter : public FilterBase,
      private util::ObjectCounter<SuperRefractionCheckImpactParameter> {
 public:
  typedef SuperRefractionCheckImpactParameterParameters Parameters_;
  static const std::string classname() {return "ufo::SuperRefractionCheckImpactParameter";}

  SuperRefractionCheckImpactParameter(ioda::ObsSpace &, const Parameters_ &,
                  std::shared_ptr<ioda::ObsDataVector<int> >,
                  std::shared_ptr<ioda::ObsDataVector<float> >);
  ~SuperRefractionCheckImpactParameter();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::superrefraction;}
  std::vector<float> calcImpactParameterModel(const std::vector<float> &,
                                          const std::vector<float> &,
                                          float, float, float) const;
  Parameters_ parameters_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_SUPERREFRACTIONCHECKIMPACTPARAMETER_H_
