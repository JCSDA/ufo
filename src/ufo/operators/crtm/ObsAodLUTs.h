/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_CRTM_OBSAODLUTS_H_
#define UFO_OPERATORS_CRTM_OBSAODLUTS_H_

#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/Fortran.h"
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
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// AodLUTs observation for UFO.
class ObsAodLUTs : public ObsOperatorBase,
                    private util::ObjectCounter<ObsAodLUTs> {
 public:
  static const std::string classname() {return "ufo::ObsAodLUTs";}

  ObsAodLUTs(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsAodLUTs();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &) const override;

// Other
  const oops::Variables & requiredVars() const override {return varin_;}

  int & toFortran() {return keyOperAodLUTs_;}
  const int & toFortran() const {return keyOperAodLUTs_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperAodLUTs_;
  const ioda::ObsSpace& odb_;
  oops::Variables varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_CRTM_OBSAODLUTS_H_
