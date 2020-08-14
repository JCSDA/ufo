# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( profile_files
      EntireSampleDataHandler.cc
      EntireSampleDataHandler.h
      ProfileChecker.cc
      ProfileChecker.h
      ProfileCheckBase.cc
      ProfileCheckBase.h
      ProfileCheckBasic.cc
      ProfileCheckBasic.h
      ProfileCheckHydrostatic.cc
      ProfileCheckHydrostatic.h
      ProfileCheckInterpolation.cc
      ProfileCheckInterpolation.h
      ProfileCheckRH.cc
      ProfileCheckRH.h
      ProfileCheckSamePDiffT.cc
      ProfileCheckSamePDiffT.h
      ProfileCheckSign.cc
      ProfileCheckSign.h
      ProfileCheckUInterp.cc
      ProfileCheckUInterp.h
      ProfileCheckUnstableLayer.cc
      ProfileCheckUnstableLayer.h
      ProfileCheckValidator.cc
      ProfileCheckValidator.h
      ProfileDataHandler.cc
      ProfileDataHandler.h
      ProfileIndices.cc
      ProfileIndices.h
      ProfileStandardLevels.cc
      ProfileStandardLevels.h
      VariableNames.h
      VariableNames.cc
)

PREPEND( _p_profile_files       "profile"       ${profile_files} )

set ( profile_src_files
      ${_p_profile_files}
      PARENT_SCOPE
)
