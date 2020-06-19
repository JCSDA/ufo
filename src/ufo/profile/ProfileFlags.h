/*
 * (C) Crown copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_PROFILEFLAGS_H_
#define UFO_PROFILE_PROFILEFLAGS_H_

#include <string>
#include <unordered_map>
#include <vector>

#include "ufo/filters/ProfileConsistencyCheckParameters.h"

#include "ufo/profile/ProfileDataBase.h"

#include "ufo/utils/Flags.h"

namespace ufo {
  class ProfileIndices;
}

namespace ufo {

  /// \brief Profile QC flags
  /// These flags are modified by the various check routines
  class ProfileFlags : public ProfileDataBase {
   public:
    ProfileFlags(ioda::ObsSpace &obsdb,
                 const ProfileConsistencyCheckParameters &options,
                 const ProfileIndices &profileIndices);

    /// Fill values of all flags for one profile
    void fillProfileValues() override;

    /// Get basic check result
    bool getBasicCheckResult() {return basicCheckResult_;}

    /// Update QC flags in total sample once they have been determined for a particular profile
    void updateFlags();

    /// Set output flags for any observation that should be rejected
    void setFlagged(const size_t nlocs, const size_t nvars,
                    std::vector<std::vector<bool>> &flagged);

    /// Reject profiles with more than a certain number of errors
    void setFinalReportFlags();

    /// Return profile tFlags
    std::vector <int> &gettFlags() {return tFlags_prof_;}

    /// Return profile RHFlags
    std::vector <int> &getRHFlags() {return RHFlags_prof_;}

    /// Return profile zFlags
    std::vector <int> &getzFlags() {return zFlags_prof_;}

    /// Return profile uFlags
    std::vector <int> &getuFlags() {return uFlags_prof_;}

    /// Return profile ReportFlags
    std::vector <int> &getReportFlags() {return ReportFlags_prof_;}

    /// Return profile tObsCorrection
    std::vector <float> &gettObsCorrection() {return tObsCorrection_prof_;}

    /// Return profile zObsCorrection
    std::vector <float> &getzObsCorrection() {return zObsCorrection_prof_;}

    /// Return entire sample counter
    std::vector <int> &getCounter(const std::string &counterName) {return Counters_[counterName];}

    /// Set basic check result
    void setBasicCheckResult(bool result) {basicCheckResult_ = result;}

    /// Increment a counter for a particular profile.
    /// Used for profile-specific error counters.
    void incrementCounter(const std::string &countname);

    /// Increment a counter for a particular profile and all subsequent ones.
    /// Used for cumulative error counters.
    void incrementCounterCumul(const std::string &countname);

    /// Write out flags to obsdb
    void writeFlags();

   private:  // functions
    /// Retrieve values of all flags for entire sample
    void retrieveAllData() override;

    /// Get counter value for all profiles from obsdb if present.
    /// Otherwise initialise to 0 everywhere.
    void initialiseCounter(const std::string &countname);

    /// Update element of one vector with a value from another vector.
    template <typename T>
      void updateValueIfPresent(const std::vector <T> &vecIn, const size_t &idxIn,
                                std::vector <T> &vecOut, const size_t &idxOut)
      {
        // Ensure neither vector is empty
        if (oops::anyVectorEmpty(vecIn, vecOut)) return;
        vecOut[idxOut] = vecIn[idxIn];
      }

   private:  // members
    //=== Entire sample values ===//

    /// Entire sample tFlags
    std::vector <int> tFlags_;

    /// Entire sample RHFlags
    std::vector <int> RHFlags_;

    /// Entire sample zFlags
    std::vector <int> zFlags_;

    /// Entire sample uFlags
    std::vector <int> uFlags_;

    /// Entire sample ReportFlags
    std::vector <int> ReportFlags_;

    /// Entire sample tObsCorrection
    std::vector <float> tObsCorrection_;

    /// Entire sample zObsCorrection
    std::vector <float> zObsCorrection_;

    //=== Individual profile values ===//

    /// Individual profile tFlags
    std::vector <int> tFlags_prof_;

    /// Individual profile RHFlags
    std::vector <int> RHFlags_prof_;

    /// Individual profile zFLags
    std::vector <int> zFlags_prof_;

    /// Individual profile uFlags
    std::vector <int> uFlags_prof_;

    /// Individual profile ReportFlags
    std::vector <int> ReportFlags_prof_;

    /// Individual profile tObsCorrection
    std::vector <float> tObsCorrection_prof_;

    /// Individual profile zObsCorrection
    std::vector <float> zObsCorrection_prof_;

    //=== Counters ===//

    /// Map of error counters for each profile in the sample
    std::unordered_map <std::string, std::vector <int>> Counters_;

    /// Number of errors in each profile
    std::vector <int> NumAnyErrors_;

    /// Number of sign changes in each profile
    std::vector <int> NumSignChange_;

    /// Number of superadiabats in each profile
    std::vector <int> NumSuperadiabat_;

    /// Number of obs with same pressure/different T or Z in each profile
    std::vector <int> NumSamePErrObs_;

    /// Number of observations with interpolation errors in each profile
    std::vector <int> NumInterpErrObs_;

    /// Number of temperature interpolation errors in each profile
    std::vector <int> NumInterpErrors_;

    /// Number of T interpolation and hydrostatic errors in each profile
    std::vector <int> NumIntHydErrors_;

    /// Number of observations with hydrostatic errors in each profile
    std::vector <int> NumHydErrObs_;

    /// Number of standard level gaps including 925 hPa in each profile
    std::vector <int> Num925Miss_;

    /// Number of standard level gaps including 100 hPa in each profile
    std::vector <int> Num100Miss_;

    /// Number of other standard level gaps in each profile
    std::vector <int> NumStdMiss_;

    /// Number of profiles with cloud errors
    std::vector <int> TotCProfs_;

    /// Number of profiles with moisture errors
    std::vector <int> TotHProfs_;

    /// Number of levels with cloud errors
    std::vector <int> TotCFlags_;

    /// Number of levels with moisture errors
    std::vector <int> TotHFlags_;

    /// Number of levels above a T threshold with moisture errors
    std::vector <int> TotLFlags_;

    /// Basic check result
    bool basicCheckResult_ = true;
  };
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILEFLAGS_H_
