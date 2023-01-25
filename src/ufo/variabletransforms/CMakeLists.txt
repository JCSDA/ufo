# (C) Copyright 2017-2018 Met Office.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( variabletransforms_files
      Cal_PressureFromHeight.h
      Cal_PressureFromHeight.cc
      Cal_HeightFromPressure.h
      Cal_HeightFromPressure.cc
      Cal_Humidity.h
      Cal_Humidity.cc
      Cal_RemapScanPosition.h
      Cal_RemapScanPosition.cc
      Cal_SatBrightnessTempFromRad.h
      Cal_SatBrightnessTempFromRad.cc
      Cal_SatRadianceFromScaledRadiance.h
      Cal_SatRadianceFromScaledRadiance.cc
      Cal_SatZenithAngleCorrection.h
      Cal_SatZenithAngleCorrection.cc
      Cal_Wind.h
      Cal_Wind.cc
      Cal_ProfileHorizontalDrift.h
      Cal_ProfileHorizontalDrift.cc 
      Cal_QNHtoQFEpressure.h
      Cal_QNHtoQFEpressure.cc
      Cal_PStar.h
      Cal_PStar.cc
      Cal_PotentialTFromT.h
      Cal_PotentialTFromT.cc
      LookupTable.h
      TransformBase.cc
      TransformBase.h
      Formulas.cc
      Formulas.h
)
if( ${gsw_FOUND} )
      set ( variabletransforms_oceanfiles
            OceanConversions/OceanDensity.cc
            OceanConversions/OceanDensity.h
            OceanConversions/OceanDepthToPressure.cc
            OceanConversions/OceanDepthToPressure.h
            OceanConversions/OceanTempToTheta.cc
            OceanConversions/OceanTempToTheta.h
            OceanConversions/OceanTempToConservativeTemp.cc
            OceanConversions/OceanTempToConservativeTemp.h
            OceanConversions/OceanPracticalSalinityToAbsoluteSalinity.cc
            OceanConversions/OceanPracticalSalinityToAbsoluteSalinity.h
      )
      
      PREPEND( _p_variabletransforms_oceanfiles       "variabletransforms"       ${variabletransforms_oceanfiles} )
endif( ${gsw_FOUND} )

PREPEND( _p_variabletransforms_files       "variabletransforms"       ${variabletransforms_files} )

set ( variabletransforms_src_files
      ${_p_variabletransforms_files}
      ${_p_variabletransforms_oceanfiles}
      PARENT_SCOPE
)
