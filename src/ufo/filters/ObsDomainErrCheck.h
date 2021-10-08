/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_OBSDOMAINERRCHECK_H_
#define UFO_FILTERS_OBSDOMAINERRCHECK_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

/// Parameters controlling the operation of the ObsDomainErrCheck filter.
class ObsDomainErrCheckParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsDomainErrCheckParameters, FilterParametersBase)

 public:
  oops::RequiredParameter<float> infltparameter{"infltparameter", this};
};

/// AMSU-A scattering check and obserr inflation

class ObsDomainErrCheck : public FilterBase,
                          private util::ObjectCounter<ObsDomainErrCheck> {
 public:
  /// The type of parameters accepted by the constructor of this filter.
  /// This typedef is used by the FilterFactory.
  typedef ObsDomainErrCheckParameters Parameters_;

  static const std::string classname() {return "ufo::ObsDomainErrCheck";}

  ObsDomainErrCheck(ioda::ObsSpace &, const Parameters_ &,
                    std::shared_ptr<ioda::ObsDataVector<int> >,
                    std::shared_ptr<ioda::ObsDataVector<float> >);
  ~ObsDomainErrCheck();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::domain;}

  Parameters_ parameters_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_OBSDOMAINERRCHECK_H_
