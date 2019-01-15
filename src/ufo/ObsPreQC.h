/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OBSPREQC_H_
#define UFO_OBSPREQC_H_

#include <ostream>

#include "eckit/config/LocalConfiguration.h"
#include "ioda/ObsSpace.h"
#include "oops/util/Printable.h"

namespace ioda {
  class ObsVector;
}

namespace ufo {
class GeoVaLs;

class ObsPreQC : public util::Printable {
 public:
  ObsPreQC(ioda::ObsSpace &, const eckit::Configuration &);
  ~ObsPreQC();

  void priorFilter(const GeoVaLs &) const {}
  void postFilter(const ioda::ObsVector &) const {}

 private:
  void print(std::ostream &) const;

  ioda::ObsSpace & obsdb_;
  const eckit::LocalConfiguration config_;
};

}  // namespace ufo

#endif  // UFO_OBSPREQC_H_
