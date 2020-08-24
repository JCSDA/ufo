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
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

/// Domain check: AMSU-A scattering check and obserr inflation
//  that obs are within domain

// Domain is defined by metadata criteria regardless of obs value.
// If obs value is required, use ObsBoundsCheck.

// The same effect can be achieved with opposite criteria through BlackList,
// the choice is a matter of convenience or which seems more natural.

class ObsDomainErrCheck : public FilterBase,
                          private util::ObjectCounter<ObsDomainErrCheck> {
 public:
  static const std::string classname() {return "ufo::ObsDomainErrCheck";}

  ObsDomainErrCheck(ioda::ObsSpace &, const eckit::Configuration &,
                    std::shared_ptr<ioda::ObsDataVector<int> >,
                    std::shared_ptr<ioda::ObsDataVector<float> >);
  ~ObsDomainErrCheck();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::domain;}

  float parameter_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_OBSDOMAINERRCHECK_H_
