/*
 * (C) Crown copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_MODELPARAMETERS_H_
#define UFO_PROFILE_MODELPARAMETERS_H_

#include <string>
#include <vector>

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"

namespace eckit {
  class Configuration;
}

namespace ufo {

  /// \brief Options related to GeoVaLs used in the profile QC code.
  class ModelParameters : public oops::Parameters {
     OOPS_CONCRETE_PARAMETERS(ModelParameters, Parameters)

   public:
    /// Height of the upper boundary of the highest model layer.
    oops::Parameter<float> zModelTop {"zModelTop", 80000.0, this};

    /// First model rho level at which there is no geographical variation in the
    /// height.
    oops::Parameter<int> firstConstantRhoLevel {"firstConstantRhoLevel", 49, this};

    /// Values of terrain-following height coordinate (eta) on theta levels.
    oops::Parameter<std::vector<float>> etaTheta {"etaTheta",
        {0.00025, 0.0006667, 0.00125, 0.002, 0.0029167, 0.004,
            0.00525, 0.0066667, 0.00825, 0.01, 0.0119167, 0.014, 0.01625, 0.0186667,
            0.02125, 0.024, 0.0269167, 0.03, 0.03325, 0.0366667, 0.04025, 0.044,
            0.0479167, 0.052, 0.05625, 0.0606667, 0.06525, 0.07, 0.0749167, 0.08,
            0.08525, 0.0906668, 0.0962505, 0.1020017, 0.1079213, 0.1140113,
            0.1202745, 0.1267154, 0.1333406, 0.1401592, 0.1471838, 0.1544313,
            0.1619238, 0.1696895, 0.1777643, 0.1861929, 0.1950307, 0.2043451,
            0.2142178, 0.2247466, 0.236048, 0.2482597, 0.2615432, 0.2760868,
            0.2921094, 0.3098631, 0.3296378, 0.3517651, 0.3766222, 0.4046373,
            0.4362943, 0.4721379, 0.5127798, 0.5589045, 0.6112759, 0.6707432,
            0.73825, 0.8148403, 0.9016668, 1},
        this};

    /// Value of terrain-following height coordinate (eta) on rho levels.
    oops::Parameter<std::vector<float>> etaRho {"etaRho",
        {0.000125, 0.0004583, 0.0009583, 0.001625, 0.0024583,
            0.0034583, 0.004625, 0.0059583, 0.0074583, 0.009125, 0.0109583,
            0.0129583, 0.015125, 0.0174583, 0.0199583, 0.022625, 0.0254583,
            0.0284583, 0.031625, 0.0349583, 0.0384583, 0.042125, 0.0459583,
            0.0499583, 0.054125, 0.0584584, 0.0629583, 0.067625, 0.0724583,
            0.0774583, 0.082625, 0.0879584, 0.0934586, 0.0991261, 0.1049615,
            0.1109663, 0.1171429, 0.123495, 0.130028, 0.1367499, 0.1436715,
            0.1508076, 0.1581776, 0.1658067, 0.1737269, 0.1819786, 0.1906118,
            0.1996879, 0.2092815, 0.2194822, 0.2303973, 0.2421538, 0.2549014,
            0.268815, 0.2840981, 0.3009862, 0.3197505, 0.3407014, 0.3641936,
            0.3906297, 0.4204658, 0.4542161, 0.4924589, 0.5358422, 0.5850902,
            0.6410096, 0.7044966, 0.7765451, 0.8582535, 0.9508334},
      this};
  };
}  // namespace ufo

#endif  // UFO_PROFILE_MODELPARAMETERS_H_

