/*
 * (C) Copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "oops/util/CompareNVectors.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "oops/util/PropertiesOfNVectors.h"

#include "ufo/profile/ModelHeightCalculator.h"

namespace ufo {
  void CalculateModelHeight(const ModelParameters &options,
                            const float orogGeoVaLs,
                            std::vector <float> &zRhoGeoVaLs,
                            std::vector <float> &zThetaGeoVaLs)
  {
    const std::vector <float> &etaThetaGeoVaLs = options.etaTheta;
    const std::vector <float> &etaRhoGeoVaLs = options.etaRho;

    if (!oops::allVectorsSameNonZeroSize(etaThetaGeoVaLs,
                                         etaRhoGeoVaLs))
      {
        oops::Log::warning() << "At least one vector is the wrong size. "
                             << "Model height calculation will not be performed." << std::endl;
        oops::Log::warning() << "Vector sizes: "
                             << oops::listOfVectorSizes(etaThetaGeoVaLs,
                                                        etaRhoGeoVaLs)
                             << std::endl;
        return;
      }

    const float missingValueFloat = util::missingValue(missingValueFloat);
    const float zModelTop = options.zModelTop;
    const int firstConstantRhoLevel = options.firstConstantRhoLevel;
    const size_t NumModLevels = etaThetaGeoVaLs.size();

    zRhoGeoVaLs.assign(NumModLevels + 1, missingValueFloat);
    zThetaGeoVaLs.assign(NumModLevels, missingValueFloat);

    // Calculate heights of rho and theta levels using terrain-following height coordinate.
    for (int jlev = 0; jlev < NumModLevels; ++jlev) {
      zRhoGeoVaLs[jlev] = etaRhoGeoVaLs[jlev] * zModelTop;
      zThetaGeoVaLs[jlev] = etaThetaGeoVaLs[jlev] * zModelTop;
    }
    // Calculate height of extra rho level.
    zRhoGeoVaLs[NumModLevels] =
      zThetaGeoVaLs[NumModLevels - 1] * 2.0 - zRhoGeoVaLs[NumModLevels - 1];
    // Quadratic method used to account for local orography.
    for (int k = 0; k < firstConstantRhoLevel; ++k) {
      zRhoGeoVaLs[k] += orogGeoVaLs *
        std::pow(1.0 - etaRhoGeoVaLs[k] / etaRhoGeoVaLs[firstConstantRhoLevel], 2.0);
      zThetaGeoVaLs[k] += orogGeoVaLs *
        std::pow(1.0 - etaThetaGeoVaLs[k] / etaRhoGeoVaLs[firstConstantRhoLevel], 2.0);
    }
  }
}  // namespace ufo
