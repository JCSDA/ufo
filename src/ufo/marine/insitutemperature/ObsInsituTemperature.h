/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_MARINE_INSITUTEMPERATURE_OBSINSITUTEMPERATURE_H_
#define UFO_MARINE_INSITUTEMPERATURE_OBSINSITUTEMPERATURE_H_

#include <memory>
#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/marine/insitutemperature/ObsInsituTemperature.interface.h"
#include "ufo/ObsOperatorBase.h"

/// Forward declarations
namespace eckit {
  class Configuration;
}

namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;

// -----------------------------------------------------------------------------
/// InsituTemperature observation operator class
class ObsInsituTemperature : public ObsOperatorBase,
                   private util::ObjectCounter<ObsInsituTemperature> {
 public:
  static const std::string classname() {return "ufo::ObsInsituTemperature";}

  ObsInsituTemperature(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsInsituTemperature();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &) const;

// Other
  const oops::Variables & variables() const {return *varin_;}

  int & toFortran() {return keyOper_;}
  const int & toFortran() const {return keyOper_;}

 private:
  void print(std::ostream &) const;
  F90hop keyOper_;
  const ioda::ObsSpace& odb_;
  std::unique_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_MARINE_INSITUTEMPERATURE_OBSINSITUTEMPERATURE_H_
