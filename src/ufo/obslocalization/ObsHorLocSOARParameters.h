/*
 * (C) Copyright 2021-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSLOCALIZATION_OBSHORLOCSOARPARAMETERS_H_
#define UFO_OBSLOCALIZATION_OBSHORLOCSOARPARAMETERS_H_

#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/obslocalization/ObsHorLocParameters.h"

namespace ufo {

/// \brief Options controlling SOAR obs localization. Inherits
/// options from general horizontal obs localization.
class ObsHorLocSOARParameters : public ObsHorLocParameters {
  OOPS_CONCRETE_PARAMETERS(ObsHorLocSOARParameters, ObsHorLocParameters)

 public:
  /// The SOAR function decay parameter
  oops::RequiredParameter<double> SOARexpDecayH{"soar horizontal decay", this};
};

}  // namespace ufo

#endif  // UFO_OBSLOCALIZATION_OBSHORLOCSOARPARAMETERS_H_
