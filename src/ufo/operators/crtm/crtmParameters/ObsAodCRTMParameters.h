/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_CRTM_CRTMPARAMETERS_OBSAODCRTMPARAMETERS_H_
#define UFO_OPERATORS_CRTM_CRTMPARAMETERS_OBSAODCRTMPARAMETERS_H_

#include "ufo/operators/crtm/crtmParameters/ObsRadianceCRTMParameters.h"

namespace ufo
{

  /// \brief Parameters controlling the CRTM AOD Forward Operator
  class ObsAodCRTMParameters : public ObsRadianceCRTMParameters
  {
    OOPS_CONCRETE_PARAMETERS(ObsAodCRTMParameters, ObsRadianceCRTMParameters)
    //
    //  Empty class for now. The yaml parameters for ObsAodCRTM are the same
    //  as for ObsRadianceCRTM at the moment. Additional parameters that are
    //  specific to the AOD operator can be added below.
    //
  };  // end class ObsAodCRTMParameters

}  // namespace ufo

#endif  // UFO_OPERATORS_CRTM_CRTMPARAMETERS_OBSAODCRTMPARAMETERS_H_
