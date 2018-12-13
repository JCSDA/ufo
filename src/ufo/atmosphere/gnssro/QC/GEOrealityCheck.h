/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_ATMOSPHERE_GNSSRO_QC_GEOREALITYCHECK_H_
#define UFO_ATMOSPHERE_GNSSRO_QC_GEOREALITYCHECK_H_

#include <ostream>
#include <string>

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "GEOrealityCheck.interface.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;

/// GEOrealityCheck:RO geophysical reality check

class GEOrealityCheck : public util::Printable,
                        private util::ObjectCounter<GEOrealityCheck> {
 public:
  static const std::string classname() {return "ufo::GEOrealityCheck";}

  GEOrealityCheck(const ioda::ObsSpace &, const eckit::Configuration &);
  ~GEOrealityCheck();
 
  void priorFilter(const GeoVaLs &) const;
  void postFilter(const ioda::ObsVector &) const;

 private:
  void print(std::ostream &) const;
  F90georealitycheck key_;
};

}  // namespace ufo

#endif  // UFO_ATMOSPHERE_GNSSRO_QC_GEOREALITYCHECK_H_
