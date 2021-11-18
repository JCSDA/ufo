/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_CRTM_CRTMPARAMETERS_OBSAODCRTMPARAMETERS_H_
#define UFO_CRTM_CRTMPARAMETERS_OBSAODCRTMPARAMETERS_H_

#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/crtm/crtmParameters/ObsRadianceCRTMParameters.h"

namespace ufo
{

  /// \brief Parameters controlling the CRTM Radiance Forward Operator
  class ObsAodCRTMParameters : public ObsRadianceCRTMParameters
  {
    ///
    ///  Empty class for now. The yaml parameters for ObsAodCRTM are
    ///  the same at the moment. Additional parameters that are
    ///  specific to the AOD operator can be added below.
    ///
  };  //  end class ObsRadianceCRTMParameters

}  // namespace ufo

#endif  // UFO_CRTM_CRTMPARAMETERS_OBSAODCRTMPARAMETERS_H_
