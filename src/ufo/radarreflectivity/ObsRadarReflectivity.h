/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_RADARREFLECTIVITY_OBSRADARREFLECTIVITY_H_
#define UFO_RADARREFLECTIVITY_OBSRADARREFLECTIVITY_H_

#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/ObsOperatorBase.h"
#include "ufo/radarreflectivity/ObsRadarReflectivity.interface.h"

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
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// RadarReflectivity observation operator class
class ObsRadarReflectivity : public ObsOperatorBase,
                   private util::ObjectCounter<ObsRadarReflectivity> {
 public:
  static const std::string classname() {return "ufo::ObsRadarReflectivity";}

  ObsRadarReflectivity(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsRadarReflectivity();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &) const override;

// Other
  const oops::Variables & variables() const {return varin_;}

  int & toFortran() {return keyOper_;}
  const int & toFortran() const {return keyOper_;}

 private:
  void print(std::ostream &) const;
  F90hop keyOper_;
  const ioda::ObsSpace& odb_;
  oops::Variables varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_RADARREFLECTIVITY_OBSRADARREFLECTIVITY_H_
