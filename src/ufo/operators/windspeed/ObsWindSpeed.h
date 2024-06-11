/*
 * (C) Crown Copyright 2024 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_WINDSPEED_OBSWINDSPEED_H_
#define UFO_OPERATORS_WINDSPEED_OBSWINDSPEED_H_

#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/ObsOperatorBase.h"
#include "ufo/ObsOperatorParametersBase.h"
#include "ufo/operators/windspeed/ObsWindSpeedParameters.h"

/// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

// -----------------------------------------------------------------------------

class WindSpeedParameters : public ObsOperatorParametersBase {
        OOPS_CONCRETE_PARAMETERS(WindSpeedParameters, ObsOperatorParametersBase)
  // NO extra parameters needed
};
/// ObsWindSpeed10m observation operator class
class ObsWindSpeed : public ObsOperatorBase,
               private util::ObjectCounter<ObsWindSpeed> {
 public:
  /// The type of parameters accepted by the constructor of this operator.
  /// This typedef is used by the ObsOperatorFactory.
  typedef ObsWindSpeedParameters Parameters_;

  static const std::string classname() {return "ufo::ObsWindSpeed";}

  // -----------------------------------------------------------------------------
  /*! \brief An operator which calculates wind speed.
  *
  * \details Calculates wind speed using background u and v wind components.
  *
  * The windspeed is calculated as:
  * \f[windspeed = \sqrt{{u}^{2}+{v}^{2}}],
  * 
  * where \f$u and \f$v are the eastward and northward wind.
  *
  * \date Mar. 2024: Created by G. Halloran (Met Office)
  *
  */
  // -----------------------------------------------------------------------------

  ObsWindSpeed(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsWindSpeed();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t&) const override;

// Other
  const oops::Variables & requiredVars() const override {return varin_;}

 private:
  void print(std::ostream &) const override;
  oops::Variable model_surface_eastward_wind_;
  oops::Variable model_surface_northward_wind_;
  oops::Variables varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_WINDSPEED_OBSWINDSPEED_H_
