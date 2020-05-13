/*
 * (C) Crown copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_PROFILECHECKVALIDATOR_H_
#define UFO_PROFILE_PROFILECHECKVALIDATOR_H_

#include <string>
#include <vector>

#include "ufo/filters/ProfileConsistencyCheckParameters.h"

#include "ufo/profile/ProfileDataBase.h"

namespace ufo {
  class ProfileIndices;
}

namespace ufo {

  /// \brief Profile QC check validator
  ///
  /// Checks various flags and intermediate values against the values produced in the OPS code.
  class ProfileCheckValidator : public ProfileDataBase {
   public:
    ProfileCheckValidator(ioda::ObsSpace &obsdb,
                          const ProfileConsistencyCheckParameters &options,
                          const ProfileIndices &profileIndices);

    /// Fill values of all variables for one profile
    void fillProfileValues() override;

    /// Validate check results against OPS values
    void validate();

    /// Get number of mismatches between values produced in this code and the OPS equivalents
    int getMismatches() const {return nMismatches_;}

    //=== Set values obtained by the profile checks ===//

    /// Set tFlags
    void settFlags(const std::vector <int> &v) {tFlags_prof_ = v;}

    /// Set zFlags
    void setzFlags(const std::vector <int> &v) {zFlags_prof_ = v;}

    /// Set ReportFlags
    void setReportFlags(const std::vector <int> &v) {ReportFlags_prof_ = v;}

    /// Set NumAnyErrors
    void setNumAnyErrors(const std::vector <int> &v) {NumAnyErrors_ = v;}

    /// Set NumSamePErrObs
    void setNumSamePErrObs(const std::vector <int> &v) {NumSamePErrObs_ = v;}

    /// Set NumSignChange
    void setNumSignChange(const std::vector <int> &v) {NumSignChange_ = v;}

    /// Set NumSuperadiabat
    void setNumSuperadiabat(const std::vector <int> &v) {NumSuperadiabat_ = v;}

    /// Set NumInterpError
    void setNumInterpErrors(const std::vector <int> &v) {NumInterpErrors_ = v;}

    /// Set NumInterpErrObs
    void setNumInterpErrObs(const std::vector <int> &v) {NumInterpErrObs_ = v;}

    /// Set Num925Miss
    void setNum925Miss(const std::vector <int> &v) {Num925Miss_ = v;}

    /// Set Num100Miss
    void setNum100Miss(const std::vector <int> &v) {Num100Miss_ = v;}

    /// Set NumStdMiss
    void setNumStdMiss(const std::vector <int> &v) {NumStdMiss_ = v;}

    /// Set NumHydErrObs
    void setNumHydErrObs(const std::vector <int> &v) {NumHydErrObs_ = v;}

    /// Set NumIntHydErrors
    void setNumIntHydErrors(const std::vector <int> &v) {NumIntHydErrors_ = v;}

    /// Set NumStd
    void setNumStd(const int &v) {NumStd_ = v;}

    /// Set NumSig
    void setNumSig(const int &v) {NumSig_ = v;}

    /// Set StdLev
    void setStdLev(const std::vector <int> &v) {StdLev_prof_ = v;}

    /// Set PBottom
    void setPBottom(const float &v) {PBottom_ = v;}

    /// Set SigAbove
    void setSigAbove(const std::vector <int> &v) {SigAbove_prof_ = v;}

    /// Set SigBelow
    void setSigBelow(const std::vector <int> &v) {SigBelow_prof_ = v;}

    /// Set IndStd
    void setIndStd(const std::vector <int> &v) {IndStd_prof_ = v;}

    /// Set LevErrors
    void setLevErrors(const std::vector <int> &v) {LevErrors_prof_ = v;}

    /// Set tInterp
    void settInterp(const std::vector <float> &v) {tInterp_prof_ = v;}

    /// Set LogP
    void setLogP(const std::vector <float> &v) {LogP_prof_ = v;}

    /// Set DC
    void setDC(const std::vector <float> &v) {DC_prof_ = v;}

    /// Set ETol
    void setETol(const std::vector <float> &v) {ETol_prof_ = v;}

    /// Set D
    void setD(const std::vector <float> &v) {D_prof_ = v;}

    /// Set E
    void setE(const std::vector <float> &v) {E_prof_ = v;}

    /// Set HydError
    void setHydError(const std::vector <int> &v) {HydError_prof_ = v;}

   private:  // functions
    /// Retrieve values of all variables for entire sample
    void retrieveAllData() override;

    /// Determine difference between two values within a certain tolerance
    template <typename T>
      bool differenceWithinTol(const T A, const T B, const float tol = 1e-10) const
      {return (std::fabs(A - B) <= tol);}

    /// Compare values with specified offset and tolerance
    template <typename T>
      void compareOutput(const std::string &desc,
                         const T val1,
                         const T val2,
                         const int offset,
                         const float tol,
                         int &n);

    /// Compare vectors of values with specified offset and tolerance
    template <typename T>
      void compareOutput(const std::string &desc,
                         const std::vector <T> &vec1,
                         const std::vector <T> &vec2,
                         const int offset,
                         const float tol,
                         int &n);

   private:  // members
    /// Number of mismatches between this code and OPS (separate for each profile)
    int nMismatches_;

    //=== Entire sample values ===//

    /// Entire sample OPS tFlags
    std::vector <int> OPS_tFlags_;

    /// Entire sample OPS zFlags
    std::vector <int> OPS_zFlags_;

    /// Entire sample OPS ReportFlags
    std::vector <int> OPS_ReportFlags_;

    /// Entire sample OPS NumStd
    std::vector <int> OPS_NumStd_;

    /// Entire sample OPS NumSig
    std::vector <int> OPS_NumSig_;

    /// Entire sample OPS NumAnyErrors
    std::vector <int> OPS_NumAnyErrors_;

    /// Entire sample OPS NumSamePErrObs
    std::vector <int> OPS_NumSamePErrObs_;

    /// Entire sample OPS NumInterpErrors
    std::vector <int> OPS_NumInterpErrors_;

    /// Entire sample OPS NumInterpErrors
    std::vector <int> OPS_NumInterpErrObs_;

    /// Entire sample OPS NumHydErrObs
    std::vector <int> OPS_NumHydErrObs_;

    /// Entire sample OPS Num925Miss
    std::vector <int> OPS_Num925Miss_;

    /// Entire sample OPS Num100Miss
    std::vector <int> OPS_Num100Miss_;

    /// Entire sample OPS NumStdMiss
    std::vector <int> OPS_NumStdMiss_;

    /// Entire sample OPS NumSignChange
    std::vector <int> OPS_NumSignChange_;

    /// Entire sample OPS NumSuperadiabat
    std::vector <int> OPS_NumSuperadiabat_;

    /// Entire sample OPS NumIntHydErrors
    std::vector <int> OPS_NumIntHydErrors_;

    /// Entire sample OPS bottom pressure
    std::vector <float> OPS_PBottom_;

    /// Entire sample OPS StdLev
    std::vector <int> OPS_StdLev_;

    /// Entire sample OPS SigBelow
    std::vector <int> OPS_SigBelow_;

    /// Entire sample OPS SigAbove
    std::vector <int> OPS_SigAbove_;

    /// Entire sample OPS LevErrors
    std::vector <int> OPS_LevErrors_;

    /// Entire sample OPS IndStd
    std::vector <int> OPS_IndStd_;

    /// Entire sample OPS tInterp
    std::vector <float> OPS_tInterp_;

    /// Entire sample OPS LogP
    std::vector <float> OPS_LogP_;

    /// Entire sample OPS DC
    std::vector <float> OPS_DC_;

    /// Entire sample OPS ETol
    std::vector <float> OPS_ETol_;

    /// Entire sample OPS D
    std::vector <float> OPS_D_;

    /// Entire sample OPS E
    std::vector <float> OPS_E_;

    /// Entire sample OPS HydError
    std::vector <int> OPS_HydError_;

    //=== Profile values ===//

    /// Individual profile OPS tFlags
    std::vector <int> OPS_tFlags_prof_;

    /// Individual profile OPS zFlags
    std::vector <int> OPS_zFlags_prof_;

    /// Individual profile OPS ReportFlags
    std::vector <int> OPS_ReportFlags_prof_;

    /// Individual profile OPS NumStd
    std::vector <int> OPS_NumStd_prof_;

    /// Individual profile OPS NumSig
    std::vector <int> OPS_NumSig_prof_;

    /// Individual profile OPS NumAnErrors
    std::vector <int> OPS_NumAnyErrors_prof_;

    /// Individual profile OPS NumSamePErrObs
    std::vector <int> OPS_NumSamePErrObs_prof_;

    /// Individual profile OPS NumInterpErrors
    std::vector <int> OPS_NumInterpErrors_prof_;

    /// Individual profile OPS NumInterpErrObs
    std::vector <int> OPS_NumInterpErrObs_prof_;

    /// Individual profile OPS NumHydErrObs
    std::vector <int> OPS_NumHydErrObs_prof_;

    /// Individual profile OPS Num925Miss
    std::vector <int> OPS_Num925Miss_prof_;

    /// Individual profile OPS Num100Miss
    std::vector <int> OPS_Num100Miss_prof_;

    /// Individual profile OPS NumStdMiss
    std::vector <int> OPS_NumStdMiss_prof_;

    /// Individual profile OPS NumSignChange
    std::vector <int> OPS_NumSignChange_prof_;

    /// Individual profile OPS NumSuperadiabat
    std::vector <int> OPS_NumSuperadiabat_prof_;

    /// Individual profile OPS NumIntHydErrors
    std::vector <int> OPS_NumIntHydErrors_prof_;

    /// Individual profile OPS bottom pressure
    std::vector <float> OPS_PBottom_prof_;

    /// Individual profile OPS StdLev
    std::vector <int> OPS_StdLev_prof_;

    /// Individual profile OPS SigBelow
    std::vector <int> OPS_SigBelow_prof_;

    /// Individual profile OPS SigAbove
    std::vector <int> OPS_SigAbove_prof_;

    /// Individual profile OPS LevErrors
    std::vector <int> OPS_LevErrors_prof_;

    /// Individual profile OPS IndStd
    std::vector <int> OPS_IndStd_prof_;

    /// Individual profile OPS tInterp
    std::vector <float> OPS_tInterp_prof_;

    /// Individual profile OPS LogP
    std::vector <float> OPS_LogP_prof_;

    /// Individual profile OPS DC
    std::vector <float> OPS_DC_prof_;

    /// Individual profile OPS ETol
    std::vector <float> OPS_ETol_prof_;

    /// Individual profile OPS D
    std::vector <float> OPS_D_prof_;

    /// Individual profile OPS E
    std::vector <float> OPS_E_prof_;

    /// Individual profile OPS HydError
    std::vector <int> OPS_HydError_prof_;

    //=== Values obtained by the profile checks ===//

    /// Individual profile tFlags
    std::vector <int> tFlags_prof_;

    /// Individual profile zFlags
    std::vector <int> zFlags_prof_;

    /// Individual profile ReportFlags
    std::vector <int> ReportFlags_prof_;

    /// Individual profile NumStd
    int NumStd_;

    /// Individual profile NumSig
    int NumSig_;

    /// Individual profile NumAnyErrors
    std::vector <int> NumAnyErrors_;

    /// Individual profile NumSamePErrObs
    std::vector <int> NumSamePErrObs_;

    /// Individual profile NumInterpErrors
    std::vector <int> NumInterpErrors_;

    /// Individual profile NumInterpErrObs
    std::vector <int> NumInterpErrObs_;

    /// Individual profile NumHydErrObs
    std::vector <int> NumHydErrObs_;

    /// Individual profile Num925Miss
    std::vector <int> Num925Miss_;

    /// Individual profile Num100Miss
    std::vector <int> Num100Miss_;

    /// Individual profile NumStdMiss
    std::vector <int> NumStdMiss_;

    /// Individual profile NumSignChange
    std::vector <int> NumSignChange_;

    /// Individual profile NumSuperadiabat
    std::vector <int> NumSuperadiabat_;

    /// Individual profile NumIntHydErrors
    std::vector <int> NumIntHydErrors_;

    /// Individual profile bottom pressure
    float PBottom_;

    /// Individual profile StdLev
    std::vector <int> StdLev_prof_;

    /// Individual profile SigBelow
    std::vector <int> SigBelow_prof_;

    /// Individual profile SigAbove
    std::vector <int> SigAbove_prof_;

    /// Individual profile LevErrors
    std::vector <int> LevErrors_prof_;

    /// Individual profile IndStd
    std::vector <int> IndStd_prof_;

    /// Individual profile tInterp
    std::vector <float> tInterp_prof_;

    /// Individual profile LogP
    std::vector <float> LogP_prof_;

    /// Individual profile DC
    std::vector <float> DC_prof_;

    /// Individual profile ETol
    std::vector <float> ETol_prof_;

    /// Individual profile D
    std::vector <float> D_prof_;

    /// Individual profile E
    std::vector <float> E_prof_;

    /// Individual profile HydError
    std::vector <int> HydError_prof_;
  };
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILECHECKVALIDATOR_H_
