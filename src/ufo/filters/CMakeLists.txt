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
      MWCLWCheck.cc
      MWCLWCheck.h
      MWSeaIceCheck.cc
      MWSeaIceCheck.h
      ObsDomainCheck.cc
      ObsDomainCheck.h
      ObsErrInflationCheck.cc
      ObsErrInflationCheck.h
      ObsDomainErrCheck.cc
      ObsDomainErrCheck.h
      ObsFilterData.cc
      ObsFilterData.h
      PreQC.cc
      PreQC.h
      QCflags.h
      QCmanager.cc
      QCmanager.h
      Thinning.cc
      Thinning.h
      Gaussian_Thinning.cc
      Gaussian_Thinning.h
      processWhere.cc
      processWhere.h
      actions/FilterAction.cc
      actions/FilterAction.h
      actions/FilterActionBase.cc
      actions/FilterActionBase.h
      actions/InflateError.cc
      actions/InflateError.h
      actions/RejectObs.cc
      actions/RejectObs.h
      obsfunctions/ObsFunction.cc
      obsfunctions/ObsFunction.h
      obsfunctions/ObsFunctionBase.cc
      obsfunctions/ObsFunctionBase.h
      obsfunctions/ObsFunctionErrInflationFactor.cc
      obsfunctions/ObsFunctionErrInflationFactor.h
      obsfunctions/ObsFunctionScattering.cc
      obsfunctions/ObsFunctionScattering.h
      obsfunctions/ObsFunctionVelocity.cc
      obsfunctions/ObsFunctionVelocity.h
)

PREPEND( _p_filters_files       "filters"       ${filters_files} )

set ( filters_src_files
      ${_p_filters_files}
      PARENT_SCOPE
)
