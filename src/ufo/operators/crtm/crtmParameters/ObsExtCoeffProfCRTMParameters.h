/*
 * (C) Copyright 2025-2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OPERATORS_CRTM_CRTMPARAMETERS_OBSEXTCOEFFPROFCRTMPARAMETERS_H_
#define UFO_OPERATORS_CRTM_CRTMPARAMETERS_OBSEXTCOEFFPROFCRTMPARAMETERS_H_

#include "ufo/operators/crtm/crtmParameters/ObsRadianceCRTMParameters.h"

namespace ufo
{

  /// \brief Parameters controlling the CRTM ExtCoeffProf Forward Operator
  class ObsExtCoeffProfCRTMParameters : public ObsRadianceCRTMParameters
  {
    OOPS_CONCRETE_PARAMETERS(ObsExtCoeffProfCRTMParameters, ObsRadianceCRTMParameters)
   public:
    /// number of lidar levels
    oops::RequiredParameter<int> nProfileLevels{"nProfileLevels", this};
  };  // end class ObsExtCoeffProfCRTMParameters

}  // namespace ufo

#endif  // UFO_OPERATORS_CRTM_CRTMPARAMETERS_OBSEXTCOEFFPROFCRTMPARAMETERS_H_
