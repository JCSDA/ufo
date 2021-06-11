/*
 * (C) Copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_MODELHEIGHTCALCULATOR_H_
#define UFO_PROFILE_MODELHEIGHTCALCULATOR_H_

#include <vector>

#include "ufo/profile/ModelParameters.h"

namespace ufo {
  /// \brief Calculate model heights on rho and theta levels.
  /// The calculation uses the terrain-following height coordinate (eta) and
  /// the local orography.
  ///
  /// \param[in] options: configuration options related to GeoVaLs.
  /// \param[in] orogGeoVaLs: orography GeoVaLs.
  /// \param[out] zRhoGeoVaLs: model heights on rho levels.
  /// \param[out] zThetaGeoVaLs: model heights on theta levels.
  void CalculateModelHeight(const ModelParameters &options,
                            const float orogGeoVaLs,
                            std::vector <float> &zRhoGeoVaLs,
                            std::vector <float> &zThetaGeoVaLs);
}  // namespace ufo

#endif  // UFO_PROFILE_MODELHEIGHTCALCULATOR_H_
