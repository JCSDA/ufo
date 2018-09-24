/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_MARINE_SEAICETHICKNESS_OBSSEAICETHICKNESSTLAD_H_
#define UFO_MARINE_SEAICETHICKNESS_OBSSEAICETHICKNESSTLAD_H_

#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/LinearObsOperatorBase.h"
#include "ufo/marine/FortranMarine.h"

// Forward declarations
namespace util {
  class DateTime;
}

namespace ioda {
  class ObsVector;
  class ObsSpace;
}

namespace ufo {
  class GeoVaLs;
  class ObsBias;
  class ObsBiasIncrement;

// -----------------------------------------------------------------------------
/// Sea-ice fraction observation for  model.
class ObsSeaIceThicknessTLAD : public LinearObsOperatorBase,
                               private util::ObjectCounter<ObsSeaIceThicknessTLAD> {
 public:
  static const std::string classname() {return "ufo::ObsSeaIceThicknessTLAD";}

  ObsSeaIceThicknessTLAD(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsSeaIceThicknessTLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, const ObsBias &);
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &, const ObsBiasIncrement &) const;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &, ObsBiasIncrement &) const;

  // Other
  const oops::Variables & variables() const {return *varin_;}

  int & toFortran() {return keyOperSeaIceThickness_;}
  const int & toFortran() const {return keyOperSeaIceThickness_;}

 private:
  void print(std::ostream &) const;
  F90hop keyOperSeaIceThickness_;
  boost::scoped_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_MARINE_SEAICETHICKNESS_OBSSEAICETHICKNESSTLAD_H_
