/*
 * (C) Crown copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_PROFILECHECKVALIDATOR_H_
#define UFO_PROFILE_PROFILECHECKVALIDATOR_H_

#include <map>
#include <set>
#include <string>
#include <vector>

#include "ufo/filters/ConventionalProfileProcessingParameters.h"

#include "ufo/profile/ProfileDataHandler.h"
#include "ufo/profile/VariableNames.h"

namespace ufo {

  /// \brief Profile QC check validator.
  ///
  /// Compares various values produced in this code to
  /// the equivalent values produced in the OPS code.
  class ProfileCheckValidator {
   public:
    explicit ProfileCheckValidator(const ConventionalProfileProcessingParameters &options);

    /// Validate check results against OPS values.
    void validate(ProfileDataHandler &profileDataHandler,
                  size_t commSize);

    /// Get number of mismatches between values produced in this code and the OPS equivalents.
    int getMismatches() const {return nMismatches_;}

   private:  // functions
    /// Compare values with specified offset and tolerance.
    template <typename T>
      void compareOutput(const std::string &desc,
                         const T val1,
                         const T val2,
                         const int offset,
                         const float tol,
                         int &n);

    /// Compare vectors of values with specified offset and tolerance.
    template <typename T>
      void compareOutput(const std::string &desc,
                         const std::vector <T> &vec1,
                         const std::vector <T> &vec2,
                         const int offset,
                         const float tol,
                         int &n);

    /// Determine difference between two values within a certain tolerance.
    template <typename T>
      bool differenceWithinTol(const T A, const T B, const float tol = 1e-10) const
      {return (std::fabs(A - B) <= tol);}

   private:  // variables
    /// Configurable parameters.
    const ConventionalProfileProcessingParameters &options_;

    /// Counters that are accumulated across profiles.
    std::map <std::string, int> cumulativeCounters_;

    /// Number of mismatches between this code and OPS (separate for each profile).
    int nMismatches_;

    /// Integer values to compare.
    std::set <std::string> valuesToCompare_int_;

    /// Float values to compare.
    std::set <std::string> valuesToCompare_float_;

    /// Values for which there is an offset between this code and OPS due to
    /// array indices starting at 0 in C++ and 1 in Fortran.
    std::map <std::string, int> comparison_offsets_;
  };
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILECHECKVALIDATOR_H_
