/*
 * (C) Copyright 2021.
 *
 * This software is developed by NOAA/NWS/EMC under the Apache 2.0 license
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 * This operator calculates modeled particulate matter (PM) at monitoring stations, such as the U.S. AirNow sites that provide PM2.5 & PM10 data.
 * With few/no code changes, it can also be applied to calculate total or/and speciated PM for applications related to other networks/datasets.
 * Unit conversion is included.
 *
 * The users are allowed to select (in the yaml file) a vertical interpolation approach to:
 * 1) match model height (above sea level = height above ground level + surface height) to station elevation - required for the AirNow application
 * or 2) match model log(air_pressure) to observation log(air_pressure) - likely suitable for use with other types of observational datasets that
 * contain pressure information, e.g. aircraft, sonde, tower..
 *  
 * Currently this tool mainly supports the calculations from the NOAA FV3-CMAQ aerosol fields (user-defined, up to 70 individual species).
 * Based on user definition in yaml, the calculation can involve the usage of the model-based scaling factors for the three modes 
 * (Aitken, Accumulation, and Coarse) that these aerosol species belong to.
 */

#ifndef UFO_INSITUPM_OBSINSITUPM_H_
#define UFO_INSITUPM_OBSINSITUPM_H_

#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"

#include "ufo/insitupm/ObsInsituPM.interface.h"
#include "ufo/insitupm/ObsInsituPMParameters.h"
#include "ufo/ObsOperatorBase.h"

/// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------
/// InsituPM observation operator class
class ObsInsituPM : public ObsOperatorBase,
                   private util::ObjectCounter<ObsInsituPM> {
 public:
  /// The type of parameters accepted by the constructor of this operator.
  /// This typedef is used by the ObsOperatorFactory.
  typedef ObsInsituPMParameters Parameters_;

  static const std::string classname() {return "ufo::ObsInsituPM";}

  ObsInsituPM(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsInsituPM();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &) const override;

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
#endif  // UFO_INSITUPM_OBSINSITUPM_H_
