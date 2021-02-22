# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.


set ( profile_files
      DataHandlerParameters.h
      EntireSampleDataHandler.cc
      EntireSampleDataHandler.h
      ModelHeightCalculator.cc
      ModelHeightCalculator.h
      ModelParameters.h
      ProfileAveragePressure.cc
      ProfileAveragePressure.h
      ProfileAverageTemperature.cc
      ProfileAverageTemperature.h
      ProfileAverageWindSpeed.cc
      ProfileAverageWindSpeed.h
      ProfileCheckBasic.cc
      ProfileChecker.cc
      ProfileChecker.h
      ProfileCheckBackgroundGeopotentialHeight.cc
      ProfileCheckBackgroundGeopotentialHeight.h
      ProfileCheckBackgroundRelativeHumidity.cc
      ProfileCheckBackgroundRelativeHumidity.h
      ProfileCheckBackgroundTemperature.cc
      ProfileCheckBackgroundTemperature.h
      ProfileCheckBackgroundWindSpeed.cc
      ProfileCheckBackgroundWindSpeed.h
      ProfileCheckBase.cc
      ProfileCheckBase.h
      ProfileCheckBasic.cc
      ProfileCheckBasic.h
      ProfileCheckHistory.cc
      ProfileCheckHistory.h
      ProfileCheckHydrostatic.cc
      ProfileCheckHydrostatic.h
      ProfileCheckInterpolation.cc
      ProfileCheckInterpolation.h
      ProfileCheckPermanentReject.cc
      ProfileCheckPermanentReject.h
      ProfileCheckRH.cc
      ProfileCheckRH.h
      ProfileCheckSamePDiffT.cc
      ProfileCheckSamePDiffT.h
      ProfileCheckSign.cc
      ProfileCheckSign.h
      ProfileCheckTime.cc
      ProfileCheckTime.h
      ProfileCheckUInterp.cc
      ProfileCheckUInterp.h
      ProfileCheckUnstableLayer.cc
      ProfileCheckUnstableLayer.h
      ProfileCheckValidator.cc
      ProfileCheckValidator.h
      ProfileDataBase.cc
      ProfileDataBase.h
      ProfileDataHandler.cc
      ProfileDataHandler.h
      ProfileDataHistory.cc
      ProfileDataHistory.h
      ProfileIndices.cc
      ProfileIndices.h
      ProfilePressure.cc
      ProfilePressure.h
      ProfileSondeFlags.cc
      ProfileSondeFlags.h
      ProfileStandardLevels.cc
      ProfileStandardLevels.h
      ProfileVerticalAveraging.cc
      ProfileVerticalAveraging.h
      ProfileVerticalInterpolation.cc
      ProfileVerticalInterpolation.h
      ProfileWindProfilerFlags.cc
      ProfileWindProfilerFlags.h
      VariableNames.cc
      VariableNames.h
)

PREPEND( _p_profile_files       "profile"       ${profile_files} )

set ( profile_src_files
      ${_p_profile_files}
      PARENT_SCOPE
)
