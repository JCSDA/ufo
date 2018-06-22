/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSINSITUTEMPERATURETLAD_H_
#define UFO_OBSINSITUTEMPERATURETLAD_H_

#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "oops/base/Variables.h"
#include "ufo/LinearObsOperatorBase.h"
#include "ioda/ObsSpace.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Logger.h"
#include "ufo/FortranMarine.h"

// Forward declarations
namespace util {
  class DateTime;
}

namespace ioda {
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsBias;
  class ObsBiasIncrement;

// -----------------------------------------------------------------------------
/// Temperature Profile observation for model.
class ObsInsituTemperatureTLAD : public LinearObsOperatorBase, 
                                 private util::ObjectCounter<ObsInsituTemperatureTLAD> {
public:
  static const std::string classname() {return "ufo::ObsInsituTemperatureTLAD";}

  ObsInsituTemperatureTLAD(const ioda::ObsSpace &, const eckit::Configuration &);    
  virtual ~ObsInsituTemperatureTLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, const ObsBias &);
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &, const ObsBiasIncrement &) const;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &, ObsBiasIncrement &) const;

  // Other
  const oops::Variables & variables() const {return *varin_;}

  int & toFortran() {return keyOperInsituTemperature_;}
  const int & toFortran() const {return keyOperInsituTemperature_;}

private:
  void print(std::ostream &) const;
  F90hop keyOperInsituTemperature_;
  const ioda::ObsSpace& odb_;
  boost::scoped_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OBSINSITUTEMPERATURETLAD_H_
