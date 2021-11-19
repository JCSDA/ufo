/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_CRTM_CRTMPARAMETERS_OBSRADIANCECRTMPARAMETERS_H_
#define UFO_CRTM_CRTMPARAMETERS_OBSRADIANCECRTMPARAMETERS_H_

#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/ObsOperatorParametersBase.h"

namespace ufo
{
  /// \brief Sub-parameters for obs options
  class CRTMObsOptionsParameters : public ObsOperatorParametersBase
  {
    OOPS_CONCRETE_PARAMETERS(CRTMObsOptionsParameters, ObsOperatorParametersBase)
   public:
    /// InspectProfileNumber
    oops::OptionalParameter<int> InspectProfile{"InspectProfileNumber", this};
    /// Sensor_ID
    oops::RequiredParameter<std::string> Sensor_ID{"Sensor_ID", this};
    /// EndianType
    oops::RequiredParameter<std::string> EndianType{"EndianType", this};
    /// CoefficientPath
    oops::RequiredParameter<std::string> CoefficientPath{"CoefficientPath", this};
    /// AerosolOption
    oops::OptionalParameter<std::string> AerosolOption{"AerosolOption", this};
    /// Surfaces
    oops::OptionalParameter<std::vector<std::string> > Surfaces{"Surfaces", this};
    /// IR Water coefficient
    oops::OptionalParameter<std::string> IRwaterCoeff{"IRwaterCoeff", "Nalli", this};
    /// VISwaterCoeff
    oops::OptionalParameter<std::string> VISwaterCoeff{"VISwaterCoeff", "NPOESS", this};
    /// IRVISlandCoeff
    oops::OptionalParameter<std::string> IRVISlandCoeff{"IRVISlandCoeff", "NPOESS", this};
    /// IRVISsnowCoeff
    oops::OptionalParameter<std::string> IRVISsnowCoeff{"IRVISsnowCoeff", "NPOESS", this};
    /// IRVISiceCoeff
    oops::OptionalParameter<std::string> IRVISiceCoeff{"IRVISiceCoeff", "NPOESS", this};
    /// MWwaterCoeff
    oops::OptionalParameter<std::string> MWwaterCoeff{"MWwaterCoeff", "FASTEM6", this};
  };  // end class ObsOptionsParameters

  /// \brief Sub-parameters for the linear obs operator
  class CRTMLinearObsOperatorParameters : public ObsOperatorParametersBase
  {
    OOPS_CONCRETE_PARAMETERS(CRTMLinearObsOperatorParameters, ObsOperatorParametersBase)
   public:
    /// Linear Absorbers
    oops::RequiredParameter<std::vector<std::string> > Absorbers{"Absorbers", this};
    /// Linear Clouds
    oops::OptionalParameter<std::vector<std::string> > Clouds{"Clouds", this};
    /// Surfaces
    oops::OptionalParameter<std::vector<std::string> > Surfaces{"Surfaces", this};
  };  // end class CRTMLinearObsOperatorParameters


  /// \brief Parameters controlling the CRTM Radiance Forward Operator
  class ObsRadianceCRTMParameters : public ObsOperatorParametersBase
  {
    OOPS_CONCRETE_PARAMETERS(ObsRadianceCRTMParameters, ObsOperatorParametersBase)
   public:
    /// Absorbers
    oops::RequiredParameter<std::vector<std::string> > Absorbers{"Absorbers", this};
    /// SurfaceWindGeoVars
    oops::Parameter<std::string> SurfaceWindGeoVars{"SurfaceWindGeoVars",
          "vector", this};
    /// Clouds
    oops::OptionalParameter<std::vector<std::string> > Clouds{"Clouds", this};
    /// Cloud_Fraction
    oops::OptionalParameter<float> Cloud_Fraction{"Cloud_Fraction", this};
    /// Salinity
    oops::OptionalParameter<bool> Salinity{"Salinity", this};
    /// Obs Options
    oops::RequiredParameter<CRTMObsOptionsParameters> obsOptions{"obs options", this};
    /// Linear Obs Operator
    oops::OptionalParameter<CRTMLinearObsOperatorParameters>
          LinearObsOperator{"linear obs operator", this};
  };  // end class ObsRadianceCRTMParameters

}  // namespace ufo

#endif  // UFO_CRTM_CRTMPARAMETERS_OBSRADIANCECRTMPARAMETERS_H_
