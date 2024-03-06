/*
 * (C) Copyright 2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_WINDUVFROMOBSSPDHOFXDIR_H_
#define UFO_FILTERS_OBSFUNCTIONS_WINDUVFROMOBSSPDHOFXDIR_H_

#include <string>

#include "oops/util/parameters/Parameter.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"
namespace ufo {

///
/// \brief An optional parameter to override the source of HofX wind components.
///
class WindUVfromObsSpdHofxDirParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(WindUVfromObsSpdHofxDirParameters, Parameters)

 public:
  /// Name of the HofX group used to replace the default group (default is HofX)
  oops::Parameter<std::string> hofxGroup{"hofx_group", "HofX", this};
};

// -----------------------------------------------------------------------------

/// \brief Compute the wind direction from the modeled velocity components and then
///        calculate the wind components based on observed wind speed and model
///        direction. For application when using observing systems that do not
///        provide window direction, such as CYGNSS
///
/// ~~~
///
/// ### Sample YAML configuration
///     - filter: Variable Assignment
///       assignments:
///       - name: DerivedObsValue/windDirection
///         type: float
///         function:
///           name: ObsFunction/WindUVfromObsSpdHofxDir
///           options:
///             hofx_group: HofX
///
class WindUVfromObsSpdHofxDir : public ObsFunctionBase<float> {
 public:
  static std::string classname() {return std::string("WindUVfromObsSpdHofxDir");}

  explicit WindUVfromObsSpdHofxDir(const eckit::LocalConfiguration &
                                 = eckit::LocalConfiguration());

  void compute(const ObsFilterData &, ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  WindUVfromObsSpdHofxDirParameters options_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_WINDUVFROMOBSSPDHOFXDIR_H_
