# (C) Copyright 2017-2018 UCAR.
# (C) Crown Copyright 2024, the Met Office.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( utils_files
      ArrowProxy.h
      Constants.h
      dataextractor/ConstrainedRange.h
      dataextractor/DataExtractor.h
      dataextractor/DataExtractor.cc
      dataextractor/DataExtractorBackend.h
      dataextractor/DataExtractorCSVBackend.h
      dataextractor/DataExtractorCSVBackend.cc
      dataextractor/DataExtractorInput.h
      dataextractor/DataExtractorNetCDFBackend.h
      dataextractor/DataExtractorNetCDFBackend.cc
      DistanceCalculator.h
      EquispacedBinSelectorBase.h
      GeodesicDistanceCalculator.h
      GeometryCalculations.cc
      GeometryCalculations.h
      IodaGroupIndices.cc
      IodaGroupIndices.h
      MaxNormDistanceCalculator.h
      metoffice/MetOfficeBMatrixStatic.cc
      metoffice/MetOfficeBMatrixStatic.h
      metoffice/MetOfficeBMatrixStatic.interface.h
      metoffice/MetOfficeBMatrixStatic.interface.F90
      metoffice/MetOfficeQCFlags.h
      metoffice/MetOfficeRMatrixRadiance.cc
      metoffice/MetOfficeRMatrixRadiance.h
      metoffice/MetOfficeRMatrixRadiance.interface.h
      metoffice/MetOfficeRMatrixRadiance.interface.F90
      metoffice/MetOfficeSort.h
      metoffice/MetOfficeObservationIDs.h
      metoffice/ufo_metoffice_bmatrixstatic_mod.f90
      metoffice/ufo_metoffice_rmatrixradiance_mod.f90
      OperatorUtils.cc
      OperatorUtils.h
      parameters/ParameterTraitsVariable.cc
      parameters/ParameterTraitsVariable.h
      PiecewiseLinearInterpolation.cc
      PiecewiseLinearInterpolation.h
      PrimitiveVariables.cc
      PrimitiveVariables.h
      ProbabilityOfGrossError.cc
      ProbabilityOfGrossError.h
      ProbabilityOfGrossErrorParameters.h
      RecordHandler.cc
      RecordHandler.h
      RecursiveSplitter.cc
      RecursiveSplitter.h
      RefractivityCalculator.F90
      RoundingEquispacedBinSelector.h
      NearestNeighborInterpolation.cc
      NearestNeighborInterpolation.h
      SaveBiasCoeffs.cc
      SaveBiasCoeffs.h
      SpatialBinSelector.h
      SpatialBinSelector.cc
      StringUtils.cc
      StringUtils.h
      SurfaceReportConstants.h
      TruncatingEquispacedBinSelector.h
      ufo_utils_mod.F90
      ufo_utils.interface.F90
      ufo_utils.interface.h
      VariableNameMap.h
      VariableNameMap.cc
      VertInterp.interface.F90
      VertInterp.interface.h
      vert_interp.F90
      thermo_utils.F90
)

if( ${gsw_FOUND} )
      set ( utils_oceanfiles
            OceanConversions/OceanConversions.interface.F90
            OceanConversions/OceanConversions.interface.h
      )

      PREPEND( _p_utils_oceanfiles       "utils"       ${utils_oceanfiles} )
endif( ${gsw_FOUND} )

PREPEND( _p_utils_files       "utils"       ${utils_files} )

set ( utils_src_files
      ${_p_utils_files}
      ${_p_utils_oceanfiles}
      PARENT_SCOPE
)
