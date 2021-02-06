/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_VARIABLETRANSFORMS_HUMIDITYRELATIVE_H_
#define UFO_FILTERS_VARIABLETRANSFORMS_HUMIDITYRELATIVE_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace formulasRH {

  // Various formulae available
  enum Method {UKMO, NOAA, Walko, Murphy, EMPTY};

  Method resolveMethods(const std::string & input);

// -------------------------------------------------------------------------------------
/*!
* \brief Calculates saturation vapor pressure over water using various methods.
*
* \b Formulae \b available:
*      - UKMO:
*        Calculation is using the Eqn 7 Sonntag (1997)
*        Reference: "Sonntag, D., Advancements in the field of hygrometry,
*        Meteorol. Zeitschrift, N. F., 3, 51-66, 1994." .
*      - Walko: 8th-order polynomial approximation (fastest choice)
*      - Murphy: most accurate
*      - NOAA: as default
*
*
* \param float temp_K: Temperature in Kelvin
* \return float saturation vapor pressure in Pascals
*/

  float SatVaporPres_fromTemp(const float temp_K,
                              const Method method = formulasRH::Method::EMPTY);
}  // namespace formulasRH

namespace ufo {

///
/// \brief One optional parameter permits overriding the default method (NOAA) used
///        to compute saturation vapor pressure.  Another optional parameter (save: true)
///        causes the calculated relative_humidity to be written back to obsdb
///
class HumidityRelativeParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(HumidityRelativeParameters, Parameters)

 public:
  /// Default: The RH calculation method can switch between default, UKMO, NOAA, NCAR,
  /// or others.  The origin or reference of each calculation can be found in source code.
  oops::Parameter<std::string> method_for_RH{"method_for_RH", "NOAA", this};
  oops::Parameter<bool> save{"save", false, this};
};

/// \brief Relative Humidity filter
///
/// The filter performs a variable conversion from specific_humidity, temperature, and
/// pressure to relative humidity. The newly calculated variable is included in the same
/// obs space. The filter can take an optional parameter to specify method for computing
/// saturation vapor pressure as a function of temperature.  Another option, save (false)
/// by default can be set true to save result of the calculation to a new variable.
///
/// Example:
///
/// \code{.yaml}
/// obs filters:
/// - filter: Relative Humidity
///   method_for_RH: UKMO             # This is NOAA by default
///   save: true                      # This is false by default
/// \endcode

class HumidityRelative : public FilterBase,
                       private util::ObjectCounter<HumidityRelative> {
 public:
  static const std::string classname() {return "ufo::HumidityRelative";}

  HumidityRelative(ioda::ObsSpace &, const eckit::Configuration &,
                  std::shared_ptr<ioda::ObsDataVector<int> >,
                  std::shared_ptr<ioda::ObsDataVector<float> >);
  ~HumidityRelative();

  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;

  const ufo::Variables & requiredVariables() const;

 private:
  void print(std::ostream &) const override;
  int qcFlag() const override {return QCflags::preQC;}
  HumidityRelativeParameters options_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_VARIABLETRANSFORMS_HUMIDITYRELATIVE_H_
