/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_CRTM_CRTMPARAMETERS_OBSRADIANCECRTMPARAMETERS_H_
#define UFO_OPERATORS_CRTM_CRTMPARAMETERS_OBSRADIANCECRTMPARAMETERS_H_

#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/GeoVaLs.h"
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
    //// use qc flags
    oops::Parameter<bool> UseQCFlagsToSkipHofX{"UseQCFlagsToSkipHofX",
       "do not calculate hofx for values not passing qc (true or false)", false, this};
    /// Sensor_ID
    oops::RequiredParameter<std::string> Sensor_ID{"Sensor_ID", this};
    /// EndianType
    oops::RequiredParameter<std::string> EndianType{"EndianType", this};
    /// CoefficientPath
    oops::RequiredParameter<std::string> CoefficientPath{"CoefficientPath", this};
    /// NetCDF CoefficientPath
    oops::OptionalParameter<std::string> NC_CoefficientPath{"NC_CoefficientPath", this};
    /// Cloud_Model
    oops::OptionalParameter<std::string> Cloud_Model{"Cloud_Model", this};
    /// CloudCoeff_File
    oops::OptionalParameter<std::string> CloudCoeff_File{"CloudCoeff_File", this};
    /// CloudCoeff_Format
    oops::OptionalParameter<std::string> CloudCoeff_Format{"CloudCoeff_Format", this};
    /// Aerosol_Model
    oops::OptionalParameter<std::string> Aerosol_Model{"Aerosol_Model", this};
    /// AerosolCoeff_File
    oops::OptionalParameter<std::string> AerosolCoeff_File{"AerosolCoeff_File", this};
    /// AerosolCoeff_Format
    oops::OptionalParameter<std::string> AerosolCoeff_Format{"AerosolCoeff_Format", this};
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
    /// Scaling factor
    oops::Parameter <double> modelUnitsCoeff{"model units coeff",
          "Conversion between model units", 1.0, this};
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
    //// use qc flags
    oops::Parameter<bool> UseQCFlagsToSkipHofX{"UseQCFlagsToSkipHofX",
       "do not calculate hofx for values not passing qc (true or false)", false, this};
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
    /// Cloud_Seeding
    oops::OptionalParameter<bool> Cloud_Seeding{"Cloud_Seeding", this};
    /// Salinity
    oops::OptionalParameter<bool> Salinity{"Salinity", this};
    /// Obs Options
    oops::RequiredParameter<CRTMObsOptionsParameters> obsOptions{"obs options", this};
    /// Linear Obs Operator
    oops::OptionalParameter<CRTMLinearObsOperatorParameters>
          LinearObsOperator{"linear obs operator", this};
    //// use qc flags
    oops::Parameter<bool> UseQCFlagsToSkipHofX{"UseQCFlagsToSkipHofX",
       "do not calculate hofx for values not passing qc (true or false)", false, this};

    /// Whether to average surface fields over sensor field of view
    oops::Parameter<bool> DoFovAverage{"do fov average", false, this};
    /// Resolution of field-of-view sampling (default arbitrary for now)
    oops::Parameter<int> FovSampleResol{"fov sample points per semimajor axis", 4, this};
    /// Write reduced GeoVaLs after FOV averaging
    oops::OptionalParameter<eckit::LocalConfiguration> gvOut{"reduced geovals output", this};
  };  // end class ObsRadianceCRTMParameters

}  // namespace ufo

#endif  // UFO_OPERATORS_CRTM_CRTMPARAMETERS_OBSRADIANCECRTMPARAMETERS_H_
