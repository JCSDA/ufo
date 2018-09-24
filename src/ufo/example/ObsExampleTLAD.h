/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

// TODO(anyone): through the file replace UFO_EXAMPLE_OBSEXAMPLETLAD_H with the unique string
#ifndef UFO_EXAMPLE_OBSEXAMPLETLAD_H_
#define UFO_EXAMPLE_OBSEXAMPLETLAD_H_

#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "ioda/ObsSpace.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/LinearObsOperatorBase.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsBias;
  class ObsBiasIncrement;

// -----------------------------------------------------------------------------
/// Example for observation operator TL and AD class
// TODO(anyone): through the file replace ObsExampleTLAD with <Your_Obs_Operator_Name>TLAD
class ObsExampleTLAD : public oops::LinearObsOperatorBase,
                       private util::ObjectCounter<ObsExampleTLAD> {
 public:
  static const std::string classname() {return "ufo::ObsExampleTLAD";}

  ObsExampleTLAD(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsExampleTLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, const ObsBias &);
  void simulateObsTL(const GeoVaLs &, ioda::ObsVector &, const ObsBiasIncrement &) const;
  void simulateObsAD(GeoVaLs &, const ioda::ObsVector &, ObsBiasIncrement &) const;

  // Other
  const oops::Variables & variables() const {return *varin_;}

  int & toFortran() {return keyOper_;}
  const int & toFortran() const {return keyOper_;}

 private:
  void print(std::ostream &) const;
  F90hop keyOper_;
  const ioda::ObsSpace& odb_;
  boost::scoped_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_EXAMPLE_OBSEXAMPLETLAD_H_
