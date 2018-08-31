/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OBSADT_H_
#define UFO_OBSADT_H_

#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "eckit/config/Configuration.h"
#include "oops/base/Variables.h"
#include "ufo/ObsOperatorBase.h"
#include "ioda/ObsSpace.h"
#include "ufo/GeoVaLs.h"
#include "ioda/Locations.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsBiasIncrement.h"
#include "ioda/ObsVector.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/marine/FortranMarine.h"

namespace ufo {

// -----------------------------------------------------------------------------
/// adt observation for UFO.
class ObsADT : public ObsOperatorBase,
               private util::ObjectCounter<ObsADT> {		    
 public:
  static const std::string classname() {return "ufo::ObsADT";}

  ObsADT(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsADT();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, const ObsBias &) const;

// Other
  const oops::Variables & variables() const {return *varin_;}

  int & toFortran() {return keyOperADT_;}
  const int & toFortran() const {return keyOperADT_;}

 private:
  void print(std::ostream &) const;
  F90hop keyOperADT_;
  const ioda::ObsSpace& odb_;
  boost::scoped_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OBSADT_H_
