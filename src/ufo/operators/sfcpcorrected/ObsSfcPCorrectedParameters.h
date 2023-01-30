/*
 * (C) Copyright 2021- UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_SFCPCORRECTED_OBSSFCPCORRECTEDPARAMETERS_H_
#define UFO_OPERATORS_SFCPCORRECTED_OBSSFCPCORRECTEDPARAMETERS_H_

#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "ufo/filters/Variable.h"
#include "ufo/ObsOperatorParametersBase.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

namespace ufo {

/// enum type for surface correction type, and ParameterTraitsHelper for it
enum class SfcPCorrectionType {
  UKMO, WRFDA, GSI
};
struct SfcPCorrectionTypeParameterTraitsHelper {
  typedef SfcPCorrectionType EnumType;
  static constexpr char enumTypeName[] = "SfcPCorrectionType";
  static constexpr util::NamedEnumerator<SfcPCorrectionType> namedValues[] = {
    { SfcPCorrectionType::UKMO, "UKMO" },
    { SfcPCorrectionType::WRFDA, "WRFDA" },
    { SfcPCorrectionType::GSI, "GSI" }
  };
};

}  // namespace ufo

namespace oops {

/// Extraction of SfcPCorrectionType parameters from config
template <>
struct ParameterTraits<ufo::SfcPCorrectionType> :
    public EnumParameterTraits<ufo::SfcPCorrectionTypeParameterTraitsHelper>
{};

}  // namespace oops

namespace ufo {

/// Configuration options recognized by the SfcPCorrected operator.
class ObsSfcPCorrectedParameters : public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsSfcPCorrectedParameters, ObsOperatorParametersBase)

 public:
  /// An optional `variables` parameter, which controls which ObsSpace
  /// variables will be simulated. This option should only be set if this operator is used as a
  /// component of the Composite operator. If `variables` is not set, the operator will simulate
  /// all ObsSpace variables. Please see the documentation of the Composite operator for further
  /// details.
  oops::OptionalParameter<std::vector<ufo::Variable>> variables{
     "variables",
     "List of variables to be simulated",
     this};

  oops::Parameter<SfcPCorrectionType> correctionType{"da_psfc_scheme",
     "Scheme used for surface pressure correction (UKMO or WRFDA)",
     SfcPCorrectionType::UKMO, this};

  /// Note: "height" default value has to be consistent with var_geomz defined
  /// in ufo_variables_mod.F90
  oops::Parameter<std::string> geovarGeomZ{"geovar_geomz",
     "Model variable for height of vertical levels",
     "height", this};

  /// Note: "surface_altitude" default value has to be consistent with var_sfc_geomz
  /// in ufo_variables_mod.F90
  oops::Parameter<std::string> geovarSfcGeomZ{"geovar_sfc_geomz",
     "Model variable for surface height",
     "surface_altitude", this};

  /// Note: "station_altitude" default value is "stationElevation"
  oops::Parameter<std::string> ObsHeightName{"station_altitude", "stationElevation", this};
};

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OPERATORS_SFCPCORRECTED_OBSSFCPCORRECTEDPARAMETERS_H_
