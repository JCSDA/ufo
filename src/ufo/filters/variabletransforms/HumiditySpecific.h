/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_VARIABLETRANSFORMS_HUMIDITYSPECIFIC_H_
#define UFO_FILTERS_VARIABLETRANSFORMS_HUMIDITYSPECIFIC_H_

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

namespace ufo {

///
/// \brief One optional parameter permits overriding the default method used to compute
///        saturation vapor pressure.  Another optional parameter (save: true)
///        causes the calculated relative_humidity to be written back to obsdb.
///
class HumiditySpecificParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(HumiditySpecificParameters, Parameters)

 public:
  /// Default: The RH calculation method can switch between default, UKMO, NOAA, NCAR,
  /// or others.  The origin or reference of each calculation can be found in source code.
  oops::Parameter<std::string> method_for_RH{"method_for_RH", "NOAA", this};
  oops::Parameter<bool> save{"save", false, this};
};

/// \brief Specific Humidity filter
///
/// The filter performs a variable conversion from relative_humidity, temperature, and
/// pressure to specific humidity. The newly calculated variable is included in the same
/// obs space. The filter can take an optional parameter to specify method for computing
/// saturation vapor pressure as a function of temperature.  Another option, save (false)
/// by default can be set true to save result of the calculation to a new variable.
///
/// Example:
///
/// \code{.yaml}
/// obs filters:
/// - filter: Specific Humidity
///   method_for_RH: UKMO               # This is NOAA by default
///   save: true                        # This is false by default
/// \endcode

class HumiditySpecific : public FilterBase,
                       private util::ObjectCounter<HumiditySpecific> {
 public:
  static const std::string classname() {return "ufo::HumiditySpecific";}

  HumiditySpecific(ioda::ObsSpace &, const eckit::Configuration &,
                  std::shared_ptr<ioda::ObsDataVector<int> >,
                  std::shared_ptr<ioda::ObsDataVector<float> >);
  ~HumiditySpecific();

  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;

  const ufo::Variables & requiredVariables() const;

 private:
  void print(std::ostream &) const override;
  int qcFlag() const override {return QCflags::preQC;}
  HumiditySpecificParameters options_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_VARIABLETRANSFORMS_HUMIDITYSPECIFIC_H_
