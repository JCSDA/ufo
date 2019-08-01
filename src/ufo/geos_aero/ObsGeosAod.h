/*
 * (C) Copyright 2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_GEOS_AERO_OBSGEOSAOD_H_
#define UFO_GEOS_AERO_OBSGEOSAOD_H_

#include <ostream>
#include <string>
#include <vector>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Logger.h"


#include "ufo/geos_aero/ObsGeosAod.interface.h"
#include "ufo/ObsOperatorBase.h"

/// Forward declarations
namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;

// -----------------------------------------------------------------------------
/// GeosAod observation operator class
class ObsGeosAod : public ObsOperatorBase,
                   private util::ObjectCounter<ObsGeosAod> {
 public:
  static const std::string classname() {return "ufo::ObsGeosAod";}

  ObsGeosAod(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsGeosAod();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &) const;

// Other
// const oops::Variables & variables() const {return *varin_;}
  const oops::Variables & variables() const {return varin_;}
 
  int & toFortran() {return keyOper_;}
  const int & toFortran() const {return keyOper_;}

 private:
  void print(std::ostream &) const;
  F90hop keyOper_;
  const ioda::ObsSpace& odb_;
//  std::unique_ptr<const oops::Variables> varin_;
  oops::Variables varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_GEOS_AERO_OBSGEOSAOD_H_
