/*
 * (C) Copyright 2026 NCAR
 * (C) Copyright 2026 University of Oklahoma / Cooperative Institute for
 *                    Severe and High-Impact Weather Research and Operations (CIWRO)
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 * Authors:
 *   Jun Park (NCAR)
 *   Yunheng Wang (CIWRO / University of Oklahoma) -- contact
 */

#ifndef UFO_OBSLOCALIZATION_OBS3DLOCPARAMETERS_H_
#define UFO_OBSLOCALIZATION_OBS3DLOCPARAMETERS_H_

#include <string>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"

#include "ufo/obslocalization/ObsHorLocParameters.h"
#include "ufo/obslocalization/ObsLocalizationParametersBase.h"
#include "ufo/utils/Constants.h"

namespace ufo {

/// \brief Options controlling 3D (combined horizontal + vertical) localization
/// using a combined normalized distance: sqrt((r_H/L_H)^2 + (r_V/L_V)^2)
class Obs3DLocParameters : public ObsLocalizationParametersBase {
  OOPS_CONCRETE_PARAMETERS(Obs3DLocParameters, ObsLocalizationParametersBase)

 public:
  /// Horizontal localization lengthscale (meters)
  oops::Parameter<double> horLengthscale{"horizontal lengthscale", 0.0, this};

  /// Vertical localization lengthscale (same units as ioda vertical coordinate)
  oops::Parameter<double> vertLengthscale{"vertical lengthscale", 0.0, this};

  /// Group in the ioda file that stores the vertical coordinate
  oops::Parameter<std::string> iodaVerticalCoordinateGroup{"ioda vertical coordinate group",
                                                           "MetaData", this};

  /// Field in the ioda file that stores the vertical coordinate
  oops::Parameter<std::string> iodaVerticalCoordinate{"ioda vertical coordinate",
                       "field in the ioda file that stores vertical coordinate", this};

  /// Apply log transformation to the vertical coordinate before computing distance
  oops::Parameter<bool> logTransform{"apply log transformation", false, this};

  /// Method for horizontal search: KD-tree or brute force
  oops::Parameter<SearchMethod> searchMethod{"search method", SearchMethod::KDTREE, this};

  /// Maximum number of obs to include
  oops::OptionalParameter<int> maxnobs{"max nobs", this};

  // Earth radius in meters (shared UFO constant, mean radius of the Earth)
  static constexpr double radius_earth = Constants::mean_earth_rad_m;
};

}  // namespace ufo

#endif  // UFO_OBSLOCALIZATION_OBS3DLOCPARAMETERS_H_
