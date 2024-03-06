/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_RADARRADIALVELOCITY_OBSRADARRADIALVELOCITY_H_
#define UFO_OPERATORS_RADARRADIALVELOCITY_OBSRADARRADIALVELOCITY_H_

#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/Parameter.h"
#include "ufo/ObsOperatorParametersBase.h"

#include "ufo/ObsOperatorBase.h"
#include "ufo/operators/radarradialvelocity/ObsRadarRadialVelocity.interface.h"

/// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
class ObsRadarRadialVelocityParameters : public ObsOperatorParametersBase  {
  OOPS_CONCRETE_PARAMETERS(ObsRadarRadialVelocityParameters, ObsOperatorParametersBase )
 public:
  /// Vertical Coordinate
  oops::Parameter<std::string> VertCoord{"VertCoord", "geometric_height", this};
};

/// RadarRadialVelocity observation operator class
class ObsRadarRadialVelocity : public ObsOperatorBase,
                   private util::ObjectCounter<ObsRadarRadialVelocity> {
 public:
  typedef ObsRadarRadialVelocityParameters Parameters_;
  typedef ioda::ObsDataVector<int> QCFlags_t;

  static const std::string classname() {return "ufo::ObsRadarRadialVelocity";}

  ObsRadarRadialVelocity(const ioda::ObsSpace &, const ObsRadarRadialVelocityParameters &);
  virtual ~ObsRadarRadialVelocity();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t &) const override;

// Other
  const oops::Variables & requiredVars() const override {return varin_;}

  int & toFortran() {return keyOper_;}
  const int & toFortran() const {return keyOper_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOper_;
  const ioda::ObsSpace& odb_;
  oops::Variables varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_RADARRADIALVELOCITY_OBSRADARRADIALVELOCITY_H_
