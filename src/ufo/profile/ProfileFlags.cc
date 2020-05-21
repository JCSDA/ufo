/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/profile/ProfileFlags.h"

namespace ufo {
  ProfileFlags::ProfileFlags(ioda::ObsSpace &obsdb,
                             const ProfileConsistencyCheckParameters &options,
                             const ProfileIndices &profileIndices)
    : ProfileDataBase(obsdb, options, profileIndices)
  {
    retrieveAllData();
  }

  void ProfileFlags::retrieveAllData()
  {
    // Retrieve QC flags from obsdb
    retrieveDataVector("tFlags", "MetaData", tFlags_);
    retrieveDataVector("zFlags", "MetaData", zFlags_);
    retrieveDataVector("ReportFlags", "MetaData", ReportFlags_);

    // Retrieve corrections from obsdb (set to zero if not present)
    retrieveDataVector("tObsCorrection", "MetaData", tObsCorrection_, true);
    retrieveDataVector("zObsCorrection", "MetaData", zObsCorrection_, true);

    // Populate map of counters with (initially empty) vectors
    Counters_["NumAnyErrors"] = NumAnyErrors_;
    Counters_["NumSamePErrObs"] = NumSamePErrObs_;
    Counters_["NumSignChange"] = NumSignChange_;
    Counters_["NumSuperadiabat"] = NumSuperadiabat_;
    Counters_["NumInterpErrors"] = NumInterpErrors_;
    Counters_["NumInterpErrObs"] = NumInterpErrObs_;
    Counters_["Num925Miss"] = Num925Miss_;
    Counters_["Num100Miss"] = Num100Miss_;
    Counters_["NumStdMiss"] = NumStdMiss_;
    Counters_["NumHydErrObs"] = NumHydErrObs_;
    Counters_["NumIntHydErrors"] = NumIntHydErrors_;

    // Initialise counters
    initialiseCounter("NumAnyErrors");
    initialiseCounter("NumSamePErrObs");
    initialiseCounter("NumSignChange");
    initialiseCounter("NumSuperadiabat");
    initialiseCounter("NumInterpErrors");
    initialiseCounter("NumInterpErrObs");
    initialiseCounter("Num925Miss");
    initialiseCounter("Num100Miss");
    initialiseCounter("NumStdMiss");
    initialiseCounter("NumHydErrObs");
    initialiseCounter("NumIntHydErrors");
  }

  void ProfileFlags::fillProfileValues()
  {
    // Fill flag information for a particular profile
    fillProfileData(tFlags_, tFlags_prof_);
    fillProfileData(zFlags_, zFlags_prof_);
    fillProfileData(ReportFlags_, ReportFlags_prof_);

    // Fill correction information
    fillProfileData(tObsCorrection_, tObsCorrection_prof_);
    fillProfileData(zObsCorrection_, zObsCorrection_prof_);
  }

  void ProfileFlags::writeFlags()
  {
    // Write flags to obsdb
    putDataVector("tFlags", "MetaData", tFlags_);
    putDataVector("zFlags", "MetaData", zFlags_);
    putDataVector("ReportFlags", "MetaData", ReportFlags_);

    // Write corrections to obsdb
    putDataVector("tObsCorrection", "MetaData", tObsCorrection_);
    putDataVector("zObsCorrection", "MetaData", zObsCorrection_);

    // Write counters to obsdb
    putDataVector("NumAnyErrors", "MetaData", Counters_["NumAnyErrors"]);
    putDataVector("NumSamePErrObs", "MetaData", Counters_["NumSamePErrObs"]);
    putDataVector("NumSignChange", "MetaData", Counters_["NumSignChange"]);
    putDataVector("NumSuperadiabat", "MetaData", Counters_["NumSuperadiabat"]);
    putDataVector("NumInterpErrors", "MetaData", Counters_["NumInterpErrors"]);
    putDataVector("NumInterpErrObs", "MetaData", Counters_["NumInterpErrObs"]);
    putDataVector("Num925Miss", "MetaData", Counters_["Num925Miss"]);
    putDataVector("Num100Miss", "MetaData", Counters_["Num100Miss"]);
    putDataVector("NumStdMiss", "MetaData", Counters_["NumStdMiss"]);
    putDataVector("NumHydErrObs", "MetaData", Counters_["NumHydErrObs"]);
    putDataVector("NumIntHydErrors", "MetaData", Counters_["NumIntHydErrors"]);
  }

  // Get counter value for all profiles from obsdb if present.
  // Otherwise initialise to 0 everywhere.
  void ProfileFlags::initialiseCounter(const std::string &countname)
  {
    retrieveCounterVector(countname, "MetaData", Counters_.at(countname));
  }

  // Increment a counter for a particular profile.
  // Used for profile-specific error counters.
  void ProfileFlags::incrementCounter(const std::string &countname)
  {
    Counters_.at(countname)[jprof_]++;
  }

  // Increment a counter for a particular profile and all subsequent ones.
  // Used for cumulative error counters.
  void ProfileFlags::incrementCounterCumul(const std::string &countname)
  {
    std::vector <int> &Counter = Counters_.at(countname);
    for (size_t j = jprof_; j < obsdb_.nrecs(); ++j) {
      Counter[j]++;
    }
  }

  // Update quantities in total sample once they have been determined for a particular profile
  void ProfileFlags::updateFlags()
  {
    int idx = 0;
    for (auto profileIndex : profileIndices_.getProfileIndices()) {
      // QC flags
      tFlags_[profileIndex] = tFlags_prof_[idx];
      zFlags_[profileIndex] = zFlags_prof_[idx];
      ReportFlags_[profileIndex] = ReportFlags_prof_[idx];

      // Correction information
      tObsCorrection_[profileIndex] = tObsCorrection_prof_[idx];
      zObsCorrection_[profileIndex] = zObsCorrection_prof_[idx];
      idx++;
    }
  }

  // Set output flags for any observation that should be rejected:
  // individual elements (FlagsElem::FinalRejectFlag) or
  // entire profile (FlagsWholeObReport::FinalRejectReport).
  void ProfileFlags::setFlagged(const size_t nlocs, const size_t nvars,
                                std::vector <std::vector <bool>> &flagged)
  {
    for (size_t jloc = 0; jloc < nlocs; ++jloc) {
      if (tFlags_[jloc] & ufo::FlagsElem::FinalRejectFlag ||
          ReportFlags_[jloc] & ufo::FlagsWholeObReport::FinalRejectReport) {
        oops::Log::debug() << "Flagging observation " << jloc << std::endl;
        // Flag all variables
        for (size_t jv = 0; jv < nvars; ++jv) {
          flagged[jv][jloc] = true;
        }
      }
    }
  }

  // Reject any profile with more than a certain number of errors
  void ProfileFlags::setFinalReportFlags()
  {
    if (Counters_.at("NumAnyErrors")[jprof_] > options_.nErrorsFail.value()) {
      oops::Log::debug() << " " << Counters_.at("NumAnyErrors")[jprof_]
                         << " errors detected, whole profile rejected" << std::endl;
      for (auto repflag : ReportFlags_prof_) {
        repflag |= ufo::FlagsWholeObReport::FinalRejectReport;
      }
    }
  }
}  // namespace ufo
