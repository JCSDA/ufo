/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OBSSEAICEFRACTION_H_
#define UFO_OBSSEAICEFRACTION_H_

#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "oops/base/Variables.h"
#include "oops/interface/ObsOperatorBase.h"
#include "ObsSpace.h"
#include "UfoTrait.h"
#include "util/ObjectCounter.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace ufo {
  class GeoVaLs;
  class Locations;
  class ObsBias;
  class ObsBiasIncrement;
  class ObsVector;

// -----------------------------------------------------------------------------
/// Total ice concentration observation for UFO.

class ObsSeaIceFraction : public oops::ObsOperatorBase<UfoTrait>,
                  private util::ObjectCounter<ObsSeaIceFraction> {
 public:
  static const std::string classname() {return "ufo::ObsSeaIceFraction";}

  ObsSeaIceFraction(const ObsSpace &, const eckit::Configuration &);
  virtual ~ObsSeaIceFraction();

// Obs Operator
  void obsEquiv(const GeoVaLs &, ObsVector &, const ObsBias &) const;

// Other
  const oops::Variables & variables() const {return *varin_;}

  int & toFortran() {return keyOperSeaIceFraction_;}
  const int & toFortran() const {return keyOperSeaIceFraction_;}

 private:
  void print(std::ostream &) const;
  F90hop keyOperSeaIceFraction_;
  const ObsSpace& odb_;
  boost::scoped_ptr<const oops::Variables> varin_;
};
// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OBSSEAICEFRACTION_H_
