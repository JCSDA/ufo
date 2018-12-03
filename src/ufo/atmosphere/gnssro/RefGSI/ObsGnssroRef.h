/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_ATMOSPHERE_GNSSRO_REFGSI_OBSGNSSROREF_H_
#define UFO_ATMOSPHERE_GNSSRO_REFGSI_OBSGNSSROREF_H_

#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/atmosphere/gnssro/RefGSI/FortranRefGSI.h"
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
  class Locations;
  class ObsBias;

// -----------------------------------------------------------------------------

/// GnssroRef observation operator
class ObsGnssroRef : public ObsOperatorBase,
                      private util::ObjectCounter<ObsGnssroRef> {
 public:
  static const std::string classname() {return "ufo::ObsGnssroRef";}

  ObsGnssroRef(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsGnssroRef();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, const ObsBias &) const;

// Other
  const oops::Variables & variables() const {return *varin_;}
  const oops::Variables & observed() const {return *varout_;}
  Locations * locateObs(const util::DateTime &, const util::DateTime &) const;

  int & toFortran() {return keyOperGnssroRef_;}
  const int & toFortran() const {return keyOperGnssroRef_;}

 private:
  void print(std::ostream &) const;
  F90hop keyOperGnssroRef_;
  const ioda::ObsSpace& odb_;
  boost::scoped_ptr<const oops::Variables> varin_;
  boost::scoped_ptr<const oops::Variables> varout_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_ATMOSPHERE_GNSSRO_REFGSI_OBSGNSSROREF_H_
