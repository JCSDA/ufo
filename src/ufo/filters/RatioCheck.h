/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_RATIOCHECK_H_
#define UFO_FILTERS_RATIOCHECK_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"
#include "ufo/filters/Variable.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

/// RatioCheck filter

class RatioCheck : public FilterBase,
                        private util::ObjectCounter<RatioCheck> {
 public:
  static const std::string classname() {return "ufo::RatioCheck";}

  RatioCheck(ioda::ObsSpace &, const eckit::Configuration &,
                  std::shared_ptr<ioda::ObsDataVector<int> >,
                  std::shared_ptr<ioda::ObsDataVector<float> >);
  ~RatioCheck();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::ratioref;}
  const Variable ref_;
  const Variable val_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_RATIOCHECK_H_
