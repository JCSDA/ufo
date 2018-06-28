/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OBSSEASURFACETEMP_H_
#define UFO_OBSSEASURFACETEMP_H_

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
#include "ufo/FortranMarine.h"

namespace ufo {

// -----------------------------------------------------------------------------

class ObsSeaSurfaceTemp : public ObsOperatorBase,
                          private util::ObjectCounter<ObsSeaSurfaceTemp> {		    
 public:
  static const std::string classname() {return "ufo::ObsSeaSurfaceTemp";}

  ObsSeaSurfaceTemp(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsSeaSurfaceTemp();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, const ObsBias &) const;

// Other
  const oops::Variables & variables() const {return *varin_;}

  int & toFortran() {return keyOperSeaSurfaceTemp_;}
  const int & toFortran() const {return keyOperSeaSurfaceTemp_;}

 private:
  void print(std::ostream &) const;
  F90hop keyOperSeaSurfaceTemp_;
  const ioda::ObsSpace& odb_;
  boost::scoped_ptr<const oops::Variables> varin_;
};

}  // namespace ufo

#endif  // UFO_OBSSEASURFACETEMP_H_

