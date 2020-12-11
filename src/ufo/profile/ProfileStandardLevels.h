/*
 * (C) Crown copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_PROFILESTANDARDLEVELS_H_
#define UFO_PROFILE_PROFILESTANDARDLEVELS_H_

#include <algorithm>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/filters/ProfileConsistencyCheckParameters.h"

#include "ufo/utils/metoffice/MetOfficeQCFlags.h"

namespace ioda {
  class ObsSpace;
}

namespace ufo {

  /// \brief Calculate standard levels
  class ProfileStandardLevels {
   public:
    explicit ProfileStandardLevels(const ProfileConsistencyCheckParameters &options);
    virtual ~ProfileStandardLevels() {}

   protected:  // functions
    /// Calculate standard levels
    void calcStdLevels(const int numProfileLevels,
                       const std::vector <float> &pressures,
                       const std::vector <float> &tObs,
                       const std::vector <int> &tFlags);

    /// Compute indices of particular standard levels for the hydrostatic check
    void findHCheckStdLevs();

    /// Calculate standard levels for U and V data
    void calcStdLevelsUV(const int numProfileLevels,
                         const std::vector <float> &pressures,
                         const std::vector <float> &uObs,
                         const std::vector <float> &vObs,
                         const std::vector <int> &uFlags);

   protected:  // variables
    /// Standard levels (hPa)
    std::vector <float> StandardLevels_;

    /// Big gaps (hPa) used in interpolation check
    std::vector <float> BigGaps_;

    /// Configurable parameters
    const ProfileConsistencyCheckParameters &optionsSL_;

    /// Number of significant levels
    int NumSig_;

    /// Number of standard levels
    int NumStd_;

    /// Index of standard levels
    std::vector <int> StdLev_;

    /// Significant level below standard level
    std::vector <int> SigBelow_;

    /// Significant level above standard level
    std::vector <int> SigAbove_;

    /// Log(Pressure) - used for vertical interpolation
    std::vector <float> LogP_;

    /// Indices of standard levels
    std::vector <int> IndStd_;

    /// Standard level index closest to 925 hPa
    int Ind925_ = 1;  // default from OPS

    /// Standard level index closest to 100 hPa
    int Ind100_ = 10;  // default from OPS
  };
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILESTANDARDLEVELS_H_
