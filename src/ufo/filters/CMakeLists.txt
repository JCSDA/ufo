# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

set ( filters_files
      BackgroundCheck.cc
      BackgroundCheck.h
      BlackList.cc
      BlackList.h
      DifferenceCheck.cc
      DifferenceCheck.h
      ObsBoundsCheck.cc
      ObsBoundsCheck.h
      ObsDomainCheck.cc
      ObsDomainCheck.h
      PreQC.cc
      PreQC.h
      QCflags.h
      QCmanager.cc
      QCmanager.h
      Thinning.cc
      Thinning.h
      processWhere.cc
      processWhere.h
)

PREPEND( _p_filters_files       "filters"       ${filters_files} )

set ( filters_src_files
      ${_p_filters_files}
      PARENT_SCOPE
)
