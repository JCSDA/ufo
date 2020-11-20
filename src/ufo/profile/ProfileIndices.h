/*
 * (C) Crown copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_PROFILE_PROFILEINDICES_H_
#define UFO_PROFILE_PROFILEINDICES_H_

#include <map>
#include <memory>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "ufo/profile/DataHandlerParameters.h"

namespace ioda {
  class ObsSpace;
}

namespace ufo {

  /// \brief Determine indices of observations making up individual profiles.
  /// The indices are computed with respect to the entire sample of observations.
  /// Used to extract relevant data and flags from the entire sample.
  ///
  /// It is important to distinguish between:
  ///   -# profile indices, which indicate the position of the profile's observations
  ///      in the entire sample of observations,
  ///   -# profile numbers, which are assigned to entire profiles
  ///      and only change when a new profile is reached.
  ///
  class ProfileIndices {
   public:
    ProfileIndices(ioda::ObsSpace &obsdb,
                   const DataHandlerParameters &options,
                   const std::vector <bool> &apply);

    /// Determine indices in entire sample for the next profile.
    void updateNextProfileIndices();

    /// Return indices for the current profile.
    const std::vector <size_t> &getProfileIndices() const {return profileIndices_;}

    /// Return number of levels to which QC checks should be applied.
    int getNumProfileLevels() const {return numProfileLevels_;}

    /// Get number of current profile.
    size_t getProfileNumCurrent() const {return profileNumCurrent_;}

   private:  // functions
    // Ensure number of profiles is consistent with quantity reported by obsdb.
    void validateTotalNumProf();

   private:  // variables
    /// Observation database.
    ioda::ObsSpace &obsdb_;

    /// Configurable parameters.
    const DataHandlerParameters &options_;

    /// Observations to apply the filter to.
    const std::vector <bool> &apply_;

    /// Profile numbers for the entire sample.
    const std::vector <size_t> profileNums_;

    /// Profile index map.
    typedef std::map<std::size_t, std::vector<std::size_t>> ProfIdxMap;

    /// Profile index map iterator.
    typedef ProfIdxMap::const_iterator ProfIdxIter;

    /// Iterator over profile indices (used for sorting).
    ProfIdxMap profidx_;

    /// Iterator pointing to current profile index (initially points to beginning).
    ProfIdxIter profidx_current_;

    /// Indices for this profile.
    std::vector <size_t> profileIndices_;

    /// Number of profile levels to which QC checks should be applied.
    int numProfileLevels_;

    /// Current profile number in the sample.
    size_t profileNumCurrent_;

    /// Next profile number to find in the sample.
    size_t profileNumToFind_;

    // Current index within the entire sample.
    size_t profIndex_;
  };
}  // namespace ufo

#endif  // UFO_PROFILE_PROFILEINDICES_H_
