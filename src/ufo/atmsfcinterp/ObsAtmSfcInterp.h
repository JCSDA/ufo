/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_ATMSFCINTERP_OBSATMSFCINTERP_H_
#define UFO_ATMSFCINTERP_OBSATMSFCINTERP_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/atmsfcinterp/ObsAtmSfcInterp.interface.h"
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
  class ObsBias;

// -----------------------------------------------------------------------------
/// AtmSfcInterp observation operator class
class ObsAtmSfcInterp : public ObsOperatorBase,
                   private util::ObjectCounter<ObsAtmSfcInterp> {
 public:
  static const std::string classname() {return "ufo::ObsAtmSfcInterp";}

  ObsAtmSfcInterp(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsAtmSfcInterp();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &) const;

// Other
  const oops::Variables & variables() const {return *varin_;}

  int & toFortran() {return keyOperAtmSfcInterp_;}
  const int & toFortran() const {return keyOperAtmSfcInterp_;}

 private:
  void print(std::ostream &) const;
  F90hop keyOperAtmSfcInterp_;
  const ioda::ObsSpace& odb_;
  std::unique_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_ATMSFCINTERP_OBSATMSFCINTERP_H_
