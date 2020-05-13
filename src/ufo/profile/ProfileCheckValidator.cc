/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <memory>
#include <string>

#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"

#include "ufo/profile/ProfileCheckValidator.h"

namespace ufo {
  ProfileCheckValidator::ProfileCheckValidator(ioda::ObsSpace &obsdb,
                                               const ProfileConsistencyCheckParameters &options,
                                               const ProfileIndices &profileIndices)
    : ProfileDataBase(obsdb, options, profileIndices)
  {
    retrieveAllData();
  }

  void ProfileCheckValidator::retrieveAllData()
  {
    // Get OPS values from obsdb
    retrieveDataVector("OPS_tFlags", "MetaData", OPS_tFlags_);
    retrieveDataVector("OPS_zFlags", "MetaData", OPS_zFlags_);
    retrieveDataVector("OPS_ReportFlags", "MetaData", OPS_ReportFlags_);
    retrieveDataVector("OPS_NumStd", "MetaData", OPS_NumStd_);
    retrieveDataVector("OPS_NumSig", "MetaData", OPS_NumSig_);
    retrieveDataVector("OPS_NumAnyErrors", "MetaData", OPS_NumAnyErrors_);
    retrieveDataVector("OPS_NumSamePErrObs", "MetaData", OPS_NumSamePErrObs_);
    retrieveDataVector("OPS_NumInterpErrors", "MetaData", OPS_NumInterpErrors_);
    retrieveDataVector("OPS_NumInterpErrObs", "MetaData", OPS_NumInterpErrObs_);
    retrieveDataVector("OPS_NumHydErrObs", "MetaData", OPS_NumHydErrObs_);
    retrieveDataVector("OPS_Num925Miss", "MetaData", OPS_Num925Miss_);
    retrieveDataVector("OPS_Num100Miss", "MetaData", OPS_Num100Miss_);
    retrieveDataVector("OPS_NumStdMiss", "MetaData", OPS_NumStdMiss_);
    retrieveDataVector("OPS_NumSignChange", "MetaData", OPS_NumSignChange_);
    retrieveDataVector("OPS_NumSuperadiabat", "MetaData", OPS_NumSuperadiabat_);
    retrieveDataVector("OPS_NumIntHydErrors", "MetaData", OPS_NumIntHydErrors_);
    retrieveDataVector("OPS_PBottom", "MetaData", OPS_PBottom_);
    retrieveDataVector("OPS_StdLev", "MetaData", OPS_StdLev_);
    retrieveDataVector("OPS_SigBelow", "MetaData", OPS_SigBelow_);
    retrieveDataVector("OPS_SigAbove", "MetaData", OPS_SigAbove_);
    retrieveDataVector("OPS_LevErrors", "MetaData", OPS_LevErrors_);
    retrieveDataVector("OPS_IndStd", "MetaData", OPS_IndStd_);
    retrieveDataVector("OPS_tInterp", "MetaData", OPS_tInterp_);
    retrieveDataVector("OPS_LogP", "MetaData", OPS_LogP_);
    retrieveDataVector("OPS_DC", "MetaData", OPS_DC_);
    retrieveDataVector("OPS_ETol", "MetaData", OPS_ETol_);
    retrieveDataVector("OPS_D", "MetaData", OPS_D_);
    retrieveDataVector("OPS_E", "MetaData", OPS_E_);
    retrieveDataVector("OPS_HydError", "MetaData", OPS_HydError_);
  }

  void ProfileCheckValidator::fillProfileValues()
  {
    // Fill OPS information for a particular profile
    fillProfileData(OPS_tFlags_, OPS_tFlags_prof_);
    fillProfileData(OPS_zFlags_, OPS_zFlags_prof_);
    fillProfileData(OPS_ReportFlags_, OPS_ReportFlags_prof_);
    fillProfileData(OPS_NumStd_, OPS_NumStd_prof_);
    fillProfileData(OPS_NumSig_, OPS_NumSig_prof_);
    fillProfileData(OPS_NumAnyErrors_, OPS_NumAnyErrors_prof_);
    fillProfileData(OPS_NumSamePErrObs_, OPS_NumSamePErrObs_prof_);
    fillProfileData(OPS_NumInterpErrors_, OPS_NumInterpErrors_prof_);
    fillProfileData(OPS_NumInterpErrObs_, OPS_NumInterpErrObs_prof_);
    fillProfileData(OPS_NumHydErrObs_, OPS_NumHydErrObs_prof_);
    fillProfileData(OPS_Num925Miss_, OPS_Num925Miss_prof_);
    fillProfileData(OPS_Num100Miss_, OPS_Num100Miss_prof_);
    fillProfileData(OPS_NumStdMiss_, OPS_NumStdMiss_prof_);
    fillProfileData(OPS_NumSignChange_, OPS_NumSignChange_prof_);
    fillProfileData(OPS_NumSuperadiabat_, OPS_NumSuperadiabat_prof_);
    fillProfileData(OPS_NumIntHydErrors_, OPS_NumIntHydErrors_prof_);
    fillProfileData(OPS_PBottom_, OPS_PBottom_prof_);
    fillProfileData(OPS_StdLev_, OPS_StdLev_prof_);
    fillProfileData(OPS_SigBelow_, OPS_SigBelow_prof_);
    fillProfileData(OPS_SigAbove_, OPS_SigAbove_prof_);
    fillProfileData(OPS_LevErrors_, OPS_LevErrors_prof_);
    fillProfileData(OPS_IndStd_, OPS_IndStd_prof_);
    fillProfileData(OPS_tInterp_, OPS_tInterp_prof_);
    fillProfileData(OPS_LogP_, OPS_LogP_prof_);
    fillProfileData(OPS_DC_, OPS_DC_prof_);
    fillProfileData(OPS_ETol_, OPS_ETol_prof_);
    fillProfileData(OPS_D_, OPS_D_prof_);
    fillProfileData(OPS_E_, OPS_E_prof_);
    fillProfileData(OPS_HydError_, OPS_HydError_prof_);
  }

  /// Comparison of single values
  template <typename T>
  void ProfileCheckValidator::compareOutput(const std::string &desc,
                                            const T val1,
                                            const T val2,
                                            const int offset,
                                            const float tol,
                                            int &n)
  {
    if (!differenceWithinTol(val1, val2 + offset, tol)) {
        oops::Log::debug() << "Mismatch for " << desc << " (OPS, this code): "
                           << val1 << ", " << val2 + offset << std::endl;
        n++;
    }
  }

  /// Comparison of vectors of values
  template <typename T>
  void ProfileCheckValidator::compareOutput(const std::string &desc,
                                            const std::vector <T> &vec1,
                                            const std::vector <T> &vec2,
                                            const int offset,
                                            const float tol,
                                            int &n)
  {
    const size_t vecsize = std::min(vec1.size(), vec2.size());

    for (size_t jvec = 0; jvec < vecsize; ++jvec) {
      if (!differenceWithinTol(vec1[jvec], vec2[jvec] + offset, tol)) {
        oops::Log::debug() << "Mismatch for " << desc << "[" << jvec << "] "
                           << "(OPS, this code): " << vec1[jvec] << ", "
                           << vec2[jvec] + offset << std::endl;
        n++;
      }
    }
  }

  void ProfileCheckValidator::validate()
  {
    nMismatches_ = 0;

    float tol = options_.Comparison_Tol.value();  // Comparison tolerance

    oops::Log::debug() << " Comparing values against OPS equivalents..." << std::endl;

    std::vector <std::string> checks = options_.Checks.value();

    // Compare OPS value with this code value
    compareOutput("tFlags",
                  OPS_tFlags_prof_,
                  tFlags_prof_,
                  0, tol, nMismatches_);
    compareOutput("zFlags",
                  OPS_zFlags_prof_,
                  zFlags_prof_,
                  0, tol, nMismatches_);
    compareOutput("ReportFlags",
                  OPS_ReportFlags_prof_,
                  ReportFlags_prof_,
                  0, tol, nMismatches_);
    compareOutput("NumAnyErrors",
                  OPS_NumAnyErrors_prof_[0],
                  NumAnyErrors_[jprof_],
                  0, tol, nMismatches_);

    for (auto check : checks) {
      if (check == "Basic") {
      } else if (check == "SamePDiffT") {
        compareOutput("NumSamePErrObs",
                      OPS_NumSamePErrObs_prof_[0],
                      NumSamePErrObs_[jprof_],
                      0, tol, nMismatches_);
      } else if (check == "Sign") {
        compareOutput("NumSignChange",
                      OPS_NumSignChange_prof_[0],
                      NumSignChange_[jprof_],
                      0, tol, nMismatches_);
      } else if (check == "UnstableLayer") {
        compareOutput("NumSuperadiabat",
                      OPS_NumSuperadiabat_prof_[0],
                      NumSuperadiabat_[jprof_],
                      0, tol, nMismatches_);
        compareOutput("PBottom",
                      OPS_PBottom_prof_[0],
                      PBottom_,
                      0, tol, nMismatches_);
      } else if (check == "Interpolation") {
        compareOutput("NumInterpErrors",
                      OPS_NumInterpErrors_prof_[0],
                      NumInterpErrors_[jprof_],
                      0, tol, nMismatches_);
        compareOutput("NumInterpErrObs",
                      OPS_NumInterpErrObs_prof_[0],
                      NumInterpErrObs_[jprof_],
                      0, tol, nMismatches_);
        compareOutput("NumStd",
                      OPS_NumStd_prof_[0],
                      NumStd_,
                      0, tol, nMismatches_);
        compareOutput("NumSig",
                      OPS_NumSig_prof_[0],
                      NumSig_,
                      0, tol, nMismatches_);
        compareOutput("StdLev",
                      OPS_StdLev_prof_,
                      StdLev_prof_,
                      1, tol, nMismatches_);
        compareOutput("SigBelow",
                      OPS_SigBelow_prof_,
                      SigBelow_prof_,
                      1, tol, nMismatches_);
        compareOutput("SigAbove",
                      OPS_SigAbove_prof_,
                      SigAbove_prof_,
                      1, tol, nMismatches_);
        compareOutput("IndStd",
                      OPS_IndStd_prof_,
                      IndStd_prof_,
                      1, tol, nMismatches_);
        compareOutput("LevErrors",
                      OPS_LevErrors_prof_,
                      LevErrors_prof_,
                      1, tol, nMismatches_);
        compareOutput("tInterp",
                      OPS_tInterp_prof_,
                      tInterp_prof_,
                      0, tol, nMismatches_);
        compareOutput("LogP",
                      OPS_LogP_prof_,
                      LogP_prof_,
                      0, tol, nMismatches_);
      } else if (check == "Hydrostatic") {
        compareOutput("Num925Miss",
                      OPS_Num925Miss_prof_[0],
                      Num925Miss_[jprof_],
                      0, tol, nMismatches_);
        compareOutput("Num100Miss",
                      OPS_Num100Miss_prof_[0],
                      Num100Miss_[jprof_],
                      0, tol, nMismatches_);
        compareOutput("NumStdMiss",
                      OPS_NumStdMiss_prof_[0],
                      NumStdMiss_[jprof_],
                      0, tol, nMismatches_);
        compareOutput("NumHydErrObs",
                      OPS_NumHydErrObs_prof_[0],
                      NumHydErrObs_[jprof_],
                      0, tol, nMismatches_);
        compareOutput("NumIntHydErrors",
                      OPS_NumIntHydErrors_prof_[0],
                      NumIntHydErrors_[jprof_],
                      0, tol, nMismatches_);
        compareOutput("DC",
                      OPS_DC_prof_,
                      DC_prof_,
                      0, tol, nMismatches_);
        compareOutput("ETol",
                      OPS_ETol_prof_,
                      ETol_prof_,
                      0, tol, nMismatches_);
        compareOutput("D",
                      OPS_D_prof_,
                      D_prof_,
                      0, tol, nMismatches_);
        compareOutput("E",
                      OPS_E_prof_,
                      E_prof_,
                      0, tol, nMismatches_);
        compareOutput("HydError",
                      OPS_HydError_prof_,
                      HydError_prof_,
                      0, tol, nMismatches_);
      }
    }

    oops::Log::debug() << " ... all comparisons done ("
                       << nMismatches_ << " mismatches)" << std::endl;
  }
}  // namespace ufo


