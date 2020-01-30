# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( utils_files
      ArrowProxy.h
      Constants.h
      DistanceCalculator.h
      EquispacedBinSelector.h
      GeodesicDistanceCalculator.h
      MaxNormDistanceCalculator.h
      parameters/ParameterTraitsVariable.h
      PiecewiseLinearInterpolation.cc
      PiecewiseLinearInterpolation.h
      RecursiveSplitter.cc
      RecursiveSplitter.h
      SpatialBinSelector.h
      SpatialBinSelector.cc
      StringUtils.cc
      StringUtils.h
      VertInterp.interface.F90
      VertInterp.interface.h
      vert_interp.F90
      thermo_utils.F90
)

PREPEND( _p_utils_files       "utils"       ${utils_files} )

set ( utils_src_files
      ${_p_utils_files}
      PARENT_SCOPE
)
