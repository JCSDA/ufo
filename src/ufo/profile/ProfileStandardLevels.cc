/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/profile/ProfileStandardLevels.h"

namespace ufo {
  ProfileStandardLevels::ProfileStandardLevels(const ProfileConsistencyCheckParameters &options)
    : optionsSL_(options)
  {
    StandardLevels_ = optionsSL_.StandardLevels.value();
    BigGaps_ = optionsSL_.BigGaps.value();
  }

  void ProfileStandardLevels::calcStdLevels(const int numProfileLevels,
                                            const std::vector <float> &pressures,
                                            const std::vector <float> &tObs,
                                            const std::vector <int> &tFlags)
  {
    oops::Log::debug() << " Finding standard levels" << std::endl;

    // Reset calculated values
    NumSig_ = 0;
    NumStd_ = 0;
    StdLev_.assign(numProfileLevels, -1);
    SigBelow_.assign(numProfileLevels, -1);
    SigAbove_.assign(numProfileLevels, -1);
    LogP_.assign(numProfileLevels, 0.0);
    IndStd_.assign(numProfileLevels, -1);

    /// Missing value (float)
    const float missingValueFloat = util::missingValue(1.0f);

    int SigPrev = -1;  // Previous significant level
    int jlevStdA = 0;  // Standard level below previous significant level
    for (int jlev = 0; jlev < numProfileLevels; ++jlev) {
      // Ignore this level if it has been flagged as rejected.
      if (tFlags[jlev] & ufo::MetOfficeQCFlags::Elem::FinalRejectFlag) continue;
      if (tObs[jlev] != missingValueFloat &&
          pressures[jlev] > optionsSL_.FS_MinP.value()) {
        LogP_[jlev] = (pressures[jlev] > 0 ? std::log(pressures[jlev]) : 0.0);
        if (tFlags[jlev] & ufo::MetOfficeQCFlags::Profile::SurfaceLevelFlag) {
          // Surface
          NumStd_++;
          StdLev_[NumStd_ - 1] = jlev;
        } else if (tFlags[jlev] & ufo::MetOfficeQCFlags::Profile::StandardLevelFlag) {
          // Standard level
          NumStd_++;
          StdLev_[NumStd_ - 1] = jlev;
          SigBelow_[NumStd_ - 1] = SigPrev;
        } else {
          NumSig_++;
          for (int jlevStd = jlevStdA; jlevStd < NumStd_; ++jlevStd) {
            SigAbove_[jlevStd] = jlev;
          }
          jlevStdA = NumStd_;
          SigPrev = jlev;
        }
      }
    }

    // Calculate IndStd_ (standard level indices)
    for (int jlevstd = 0; jlevstd < NumStd_; ++jlevstd) {
      int jlev = StdLev_[jlevstd];  // Standard level
      if (tFlags[jlev] & ufo::MetOfficeQCFlags::Profile::SurfaceLevelFlag) continue;
      int IPStd = std::round(pressures[jlev] * 0.01);  // Pressure rounded to nearest hPa
      for (size_t i = 0; i < StandardLevels_.size(); ++i) {
        if (IPStd == StandardLevels_[i])
          IndStd_[jlevstd] = static_cast<int> (i);  // Index at which standard level appears
        if (StandardLevels_[i] <= IPStd) break;
      }
    }
  }

  void ProfileStandardLevels::findHCheckStdLevs()
  {
    // Find indices in StandardLevels that are closest to 925 and 100 hPa
    // Required for hydrostatic check

    // StandardLevels_ must contain decreasing pressures
    for (size_t jstdlev = 0; jstdlev < StandardLevels_.size() - 1; ++jstdlev) {
      if (StandardLevels_[jstdlev] < StandardLevels_[jstdlev + 1]) {
          throw eckit::BadValue("Standard levels in wrong order", Here());
        }
    }

    // Find indices using reverse iterators
    auto it925 = std::lower_bound(StandardLevels_.rbegin(), StandardLevels_.rend(), 925.0);
    Ind925_ = std::distance(StandardLevels_.begin(), it925.base()) - 1;
    auto it100 = std::lower_bound(StandardLevels_.rbegin(), StandardLevels_.rend(), 100.0);
    Ind100_ = std::distance(StandardLevels_.begin(), it100.base()) - 1;
  }

  void ProfileStandardLevels::calcStdLevelsUV(const int numProfileLevels,
                                              const std::vector <float> &pressures,
                                              const std::vector <float> &uObs,
                                              const std::vector <float> &vObs,
                                              const std::vector <int> &uFlags)
  {
    oops::Log::debug() << " Finding standard levels for U and V data" << std::endl;

    // Reset calculated values
    NumSig_ = 0;
    NumStd_ = 0;
    StdLev_.assign(numProfileLevels, -1);
    SigBelow_.assign(numProfileLevels, -1);
    SigAbove_.assign(numProfileLevels, -1);
    LogP_.assign(numProfileLevels, 0.0);

    /// Missing value (float)
    const float missingValueFloat = util::missingValue(1.0f);

    int SigPrev = -1;  // Previous significant level
    int jlevStdA = 0;  // Standard level below previous significant level

    for (int jlev = 0; jlev < numProfileLevels; ++jlev) {
      if (uObs[jlev] != missingValueFloat && vObs[jlev] != missingValueFloat) {
        LogP_[jlev] = std::log(pressures[jlev]);
        if (uFlags[jlev] & ufo::MetOfficeQCFlags::Profile::StandardLevelFlag) {  // Standard level
          NumStd_++;
          StdLev_[NumStd_ - 1] = jlev;
          SigBelow_[NumStd_ - 1] = SigPrev;
        } else {
          NumSig_++;
          for (int jlevstd = jlevStdA; jlevstd < NumStd_; ++jlevstd) {
            SigAbove_[jlevstd] = jlev;
          }
          jlevStdA = NumStd_;
          SigPrev = jlev;
        }
      }
    }
  }
}  // namespace ufo
