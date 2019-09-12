/*
 * (C) Copyright 2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_THINNING_H_
#define UFO_FILTERS_THINNING_H_

#include <ostream>
#include <string>
#include <vector>

#include "boost/shared_ptr.hpp"

#include "ioda/ObsDataVector.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/filters/FilterBase.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

/// Thinning: randonly thin a given percentage of observations

class Thinning : public FilterBase,
                 private util::ObjectCounter<Thinning> {
 public:
  static const std::string classname() {return "ufo::Thinning";}

  Thinning(ioda::ObsSpace &, const eckit::Configuration &,
           boost::shared_ptr<ioda::ObsDataVector<int> >,
           boost::shared_ptr<ioda::ObsDataVector<float> >);
  ~Thinning();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, std::vector<std::vector<bool>> &) const override;
};

}  // namespace ufo

#endif  // UFO_FILTERS_THINNING_H_
