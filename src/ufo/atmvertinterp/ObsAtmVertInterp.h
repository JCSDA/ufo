/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_ATMVERTINTERP_OBSATMVERTINTERP_H_
#define UFO_ATMVERTINTERP_OBSATMVERTINTERP_H_

#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/atmvertinterp/ObsAtmVertInterp.interface.h"
#include "ufo/ObsOperatorBase.h"

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

/// AtmVertInterp observation operator
class ObsAtmVertInterp : public ObsOperatorBase,
                      private util::ObjectCounter<ObsAtmVertInterp> {
 public:
  static const std::string classname() {return "ufo::ObsAtmVertInterp";}

  ObsAtmVertInterp(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsAtmVertInterp();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &) const;

// Other
  const oops::Variables & variables() const {return varin_;}

  int & toFortran() {return keyOperAtmVertInterp_;}
  const int & toFortran() const {return keyOperAtmVertInterp_;}

 private:
  void print(std::ostream &) const;
  F90hop keyOperAtmVertInterp_;
  const ioda::ObsSpace& odb_;
  oops::Variables varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_ATMVERTINTERP_OBSATMVERTINTERP_H_
