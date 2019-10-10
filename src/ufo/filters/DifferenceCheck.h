/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_DIFFERENCECHECK_H_
#define UFO_FILTERS_DIFFERENCECHECK_H_

#include <ostream>
#include <string>
#include <vector>

#include "boost/shared_ptr.hpp"

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

/// DifferenceCheck filter

class DifferenceCheck : public FilterBase,
                        private util::ObjectCounter<DifferenceCheck> {
 public:
  static const std::string classname() {return "ufo::DifferenceCheck";}

  DifferenceCheck(ioda::ObsSpace &, const eckit::Configuration &,
                  boost::shared_ptr<ioda::ObsDataVector<int> >,
                  boost::shared_ptr<ioda::ObsDataVector<float> >);
  ~DifferenceCheck();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, std::vector<std::vector<bool>> &) const override;

  const std::string ref_;
  const std::string val_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_DIFFERENCECHECK_H_
