/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_PREQC_H_
#define UFO_FILTERS_PREQC_H_

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

class PreQCParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(PreQCParameters, FilterParametersBase)

 public:
  /// The ObsSpace group holding PreQC flags. By default, 'PreQC'.
  oops::Parameter<std::string> inputQC{"inputQC", "PreQC", this};
  /// Minimum PreQC flag denoting "pass". By default, zero.
  oops::Parameter<int> minvalue{"minvalue", 0, this};
  /// Maximum PreQC flag denoting "pass". By default, zero.
  oops::Parameter<int> maxvalue{"maxvalue", 0, this};
};

class PreQC : public FilterBase, private util::ObjectCounter<PreQC> {
 public:
  /// The type of parameters accepted by the constructor of this filter.
  /// This typedef is used by the FilterFactory.
  typedef PreQCParameters Parameters_;

  static const std::string classname() {return "ufo::PreQC";}

  PreQC(ioda::ObsSpace &, const Parameters_ &,
        std::shared_ptr<ioda::ObsDataVector<int> >,
        std::shared_ptr<ioda::ObsDataVector<float> >);

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::preQC;}

  Parameters_ parameters_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_PREQC_H_
