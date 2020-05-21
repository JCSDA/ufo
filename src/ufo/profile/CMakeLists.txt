# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( profile_files
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
      ProfileCheckSamePDiffT.cc
      ProfileCheckSamePDiffT.h
      ProfileCheckSign.cc
      ProfileCheckSign.h
      ProfileCheckUnstableLayer.cc
      ProfileCheckUnstableLayer.h
      ProfileCheckValidator.cc
      ProfileCheckValidator.h
      ProfileData.cc
      ProfileData.h
      ProfileDataBase.cc
      ProfileDataBase.h
      ProfileFlags.cc
      ProfileFlags.h
      ProfileIndices.cc
      ProfileIndices.h
      ProfileStandardLevels.cc
      ProfileStandardLevels.h
)

PREPEND( _p_profile_files       "profile"       ${profile_files} )

set ( profile_src_files
      ${_p_profile_files}
      PARENT_SCOPE
)
