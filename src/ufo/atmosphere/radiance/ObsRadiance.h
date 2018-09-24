/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_ATMOSPHERE_RADIANCE_OBSRADIANCE_H_
#define UFO_ATMOSPHERE_RADIANCE_OBSRADIANCE_H_

#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/atmosphere/FortranAtmosphere.h"
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
/// Radiance observation for UFO.
class ObsRadiance : public ObsOperatorBase,
                    private util::ObjectCounter<ObsRadiance> {
 public:
  static const std::string classname() {return "ufo::ObsRadiance";}

  ObsRadiance(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsRadiance();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, const ObsBias &) const;

// Other
  const oops::Variables & variables() const {return *varin_;}

  int & toFortran() {return keyOperRadiance_;}
  const int & toFortran() const {return keyOperRadiance_;}

 private:
  void print(std::ostream &) const;
  F90hop keyOperRadiance_;
  const ioda::ObsSpace& odb_;
  boost::scoped_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_ATMOSPHERE_RADIANCE_OBSRADIANCE_H_
