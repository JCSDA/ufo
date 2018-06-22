/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OBSAOD_H_
#define UFO_OBSAOD_H_

#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/FortranConstituents.h"
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
class ObsBiasIncrement;

// -----------------------------------------------------------------------------
/// Aod observation for UFO.
class ObsAod : public ObsOperatorBase,
               private util::ObjectCounter<ObsAod> {
 public:
  static const std::string classname() {return "ufo::ObsAod";}

  ObsAod(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsAod();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, const ObsBias &) const;

// Other
  const oops::Variables & variables() const {return *varin_;}

  int & toFortran() {return keyOperAod_;}
  const int & toFortran() const {return keyOperAod_;}

 private:
  void print(std::ostream &) const;
  F90hop keyOperAod_;
  const ioda::ObsSpace& odb_;
  boost::scoped_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OBSAOD_H_
