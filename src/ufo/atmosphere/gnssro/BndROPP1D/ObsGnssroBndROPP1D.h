/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_ATMOSPHERE_GNSSRO_BNDROPP1D_OBSGNSSROBNDROPP1D_H_
#define UFO_ATMOSPHERE_GNSSRO_BNDROPP1D_OBSGNSSROBNDROPP1D_H_

#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/atmosphere/gnssro/BndROPP1D/FortranBndROPP1D.h"
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

/// GnssroBndROPP1D observation operator
class ObsGnssroBndROPP1D : public ObsOperatorBase,
                        private util::ObjectCounter<ObsGnssroBndROPP1D> {
 public:
  static const std::string classname() {return "ufo::ObsGnssroBndROPP1D";}

  ObsGnssroBndROPP1D(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsGnssroBndROPP1D();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, const ObsBias &) const;

// Other
  const oops::Variables & variables() const {return *varin_;}
  const oops::Variables & observed() const {return *varout_;}

  int & toFortran() {return keyOperGnssroBndROPP1D_;}
  const int & toFortran() const {return keyOperGnssroBndROPP1D_;}

 private:
  void print(std::ostream &) const;
  F90hop keyOperGnssroBndROPP1D_;
  const ioda::ObsSpace& odb_;
  boost::scoped_ptr<const oops::Variables> varin_;
  boost::scoped_ptr<const oops::Variables> varout_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_ATMOSPHERE_GNSSRO_BNDROPP1D_OBSGNSSROBNDROPP1D_H_
