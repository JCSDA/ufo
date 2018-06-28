/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSADTTLAD_H_
#define UFO_OBSADTTLAD_H_

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
/// ADT observation for  model.
class ObsADTTLAD : public LinearObsOperatorBase, 
                   private util::ObjectCounter<ObsADTTLAD> {
				
public:
  static const std::string classname() {return "ufo::ObsADTTLAD";}

  ObsADTTLAD(const ioda::ObsSpace &, const eckit::Configuration &);    
  virtual ~ObsADTTLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, const ObsBias &);
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &, const ObsBiasIncrement &) const;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &, ObsBiasIncrement &) const;

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
#endif  // UFO_OBSADTTLAD_H_
