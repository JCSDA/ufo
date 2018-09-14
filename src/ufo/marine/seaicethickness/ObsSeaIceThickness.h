/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_MARINE_SEAICETHICKNESS_OBSSEAICETHICKNESS_H_
#define UFO_MARINE_SEAICETHICKNESS_OBSSEAICETHICKNESS_H_

#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/marine/FortranMarine.h"
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
  class ObsBias;

// -----------------------------------------------------------------------------
/// Total ice concentration observation for UFO.
class ObsSeaIceThickness : public ObsOperatorBase,
                           private util::ObjectCounter<ObsSeaIceThickness> {
 public:
  static const std::string classname() {return "ufo::ObsSeaIceThickness";}

  ObsSeaIceThickness(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsSeaIceThickness();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, const ObsBias &) const;

// Other
  const oops::Variables & variables() const {return *varin_;}

  int & toFortran() {return keyOperSeaIceThickness_;}
  const int & toFortran() const {return keyOperSeaIceThickness_;}

 private:
  void print(std::ostream &) const;
  F90hop keyOperSeaIceThickness_;
  const ioda::ObsSpace& odb_;
  boost::scoped_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_MARINE_SEAICETHICKNESS_OBSSEAICETHICKNESS_H_
