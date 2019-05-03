/*
 * (C) Copyright 2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_THINNING_H_
#define UFO_THINNING_H_

#include <ostream>
#include <string>

#include "boost/shared_ptr.hpp"

#include "eckit/config/LocalConfiguration.h"
#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace ioda {
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;

/// Thinning: randonly thin a given percentage of observations

class Thinning : public util::Printable,
                  private util::ObjectCounter<Thinning> {
 public:
  static const std::string classname() {return "ufo::Thinning";}

  Thinning(ioda::ObsSpace &, const eckit::Configuration &,
           boost::shared_ptr<ioda::ObsDataVector<int> >,
           boost::shared_ptr<ioda::ObsDataVector<float> >);
  ~Thinning();

  void priorFilter(const GeoVaLs &) const;
  void postFilter(const ioda::ObsVector &) const {}

  const oops::Variables & requiredGeoVaLs() const {return geovars_;}

 private:
  void print(std::ostream &) const;

  ioda::ObsSpace & obsdb_;
  const eckit::LocalConfiguration config_;
  oops::Variables geovars_;
  ioda::ObsDataVector<int> & flags_;
};

}  // namespace ufo

#endif  // UFO_THINNING_H_
