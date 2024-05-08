# (C) Copyright 2020 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( predictor_files
  PredictorBase.h
  PredictorBase.cc
  CloudLiquidWater.h
  CloudLiquidWater.cc
  Constant.h
  Constant.cc
  CosineOfLatitudeTimesOrbitNode.h
  CosineOfLatitudeTimesOrbitNode.cc
  Emissivity.h
  Emissivity.cc
  InterpolateDataFromFile.h
  InterpolateDataFromFile.cc
  LapseRate.h
  LapseRate.cc
  Legendre.h
  Legendre.cc
  OrbitalAngle.h
  OrbitalAngle.cc
  ReadBias.h
  ReadBias.cc
  SatelliteSelector.h
  SatelliteSelector.cc
  ScanAngle.h
  ScanAngle.cc
  SineOfLatitude.h
  SineOfLatitude.cc
  Thickness.h
  Thickness.cc
  ObsMetaDataPredictor.h
  ObsMetaDataPredictor.cc
)

PREPEND( _p_predictor_files     "predictors"            ${predictor_files} )

set ( predictor_src_files
      ${_p_predictor_files}
      PARENT_SCOPE
)
