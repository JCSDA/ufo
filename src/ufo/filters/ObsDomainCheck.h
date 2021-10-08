/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_OBSDOMAINCHECK_H_
#define UFO_FILTERS_OBSDOMAINCHECK_H_

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

class ObsDomainCheckParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsDomainCheckParameters, FilterParametersBase)

  // This filter doesn't take any extra parameters.
};

/// Domain check: generic check that obs are within domain

// Domain is defined by metadata criteria regardless of obs value.
// If obs value is required, use ObsBoundsCheck.

// The same effect can be achieved with opposite criteria through BlackList,
// the choice is a matter of convenience or which seems more natural.

class ObsDomainCheck : public FilterBase,
                       private util::ObjectCounter<ObsDomainCheck> {
 public:
  /// The type of parameters accepted by the constructor of this filter.
  /// This typedef is used by the FilterFactory.
  typedef ObsDomainCheckParameters Parameters_;

  static const std::string classname() {return "ufo::ObsDomainCheck";}

  ObsDomainCheck(ioda::ObsSpace &, const Parameters_ &,
                 std::shared_ptr<ioda::ObsDataVector<int> >,
                 std::shared_ptr<ioda::ObsDataVector<float> >);
  ~ObsDomainCheck();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::domain;}

  Parameters_ parameters_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_OBSDOMAINCHECK_H_
