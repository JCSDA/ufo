/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_OBSBOUNDSCHECK_H_
#define UFO_FILTERS_OBSBOUNDSCHECK_H_

#include <ostream>
#include <string>
#include <vector>

#include "boost/shared_ptr.hpp"

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

/// ObsBoundsCheck: generic quality control based on observation data only

// Check that observations are within some bounds over some domain

class ObsBoundsCheck : public FilterBase,
                       private util::ObjectCounter<ObsBoundsCheck> {
 public:
  static const std::string classname() {return "ufo::ObsBoundsCheck";}

  ObsBoundsCheck(ioda::ObsSpace &, const eckit::Configuration &,
                 boost::shared_ptr<ioda::ObsDataVector<int> >,
                 boost::shared_ptr<ioda::ObsDataVector<float> >);
  ~ObsBoundsCheck();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::bounds;}
};

}  // namespace ufo

#endif  // UFO_FILTERS_OBSBOUNDSCHECK_H_
