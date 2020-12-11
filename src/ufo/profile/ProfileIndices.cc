/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <utility>

#include "ufo/profile/ProfileIndices.h"

namespace ufo {
  ProfileIndices::ProfileIndices(ioda::ObsSpace &obsdb,
                                 const DataHandlerParameters &options,
                                 const std::vector <bool> &apply)
    : obsdb_(obsdb),
      options_(options),
      apply_(apply),
      profileNums_(obsdb.recnum()),
      profileNumCurrent_(0),
      profileNumToFind_(0),
      profIndex_(0)
  {
    // If not sorting observations, ensure number of profiles is consistent
    // with quantity reported by obsdb.
    // (If sorting is imposed, nothing is assumed about the ordering of the input data
    // so this validation should not be performed.)
    if (obsdb_.obs_sort_var().empty() && options_.ValidateTotalNumProf.value()) {
      validateTotalNumProf();
    }

    // If sorting observations, point to beginning of record index iterator
    if (!obsdb_.obs_sort_var().empty() &&
        obsdb_.obs_sort_order() == "descending") {
      profidx_current_ = obsdb_.recidx_begin();
    }
  }

  void ProfileIndices::updateNextProfileIndices()
  {
    profileIndices_.clear();

    // If there are no profiles in the sample, warn and exit
    if (profileNums_.size() == 0) {
      oops::Log::debug() << "No profiles found in the sample" << std::endl;
      return;
    }

    // Determine indices in the full sample that correspond to
    // the profile ID to be found.
    // The method used to fill the profile indices depends upon the sorting chosen.
    if (obsdb_.obs_sort_var().empty()) {
      // If no sorting has been specified just increment indices
      while (profIndex_ < profileNums_.size() && profileNums_[profIndex_] == profileNumToFind_) {
        if (apply_[profIndex_]) {
          profileIndices_.emplace_back(profIndex_);
        }
        profIndex_++;
      }
    } else {
      if (obsdb_.obs_sort_order() == "descending") {
        // Sort variable (usually pressure) in descending order
        // Sorted indices for the current profile
        const std::vector<std::size_t> profidx_sorted = profidx_current_->second;
        auto it_profidx_sorted = profidx_sorted.begin();
        while (profIndex_ < profileNums_.size() && profileNums_[profIndex_] == profileNumToFind_) {
          if (apply_[profIndex_]) {
            profileIndices_.emplace_back(*it_profidx_sorted);
          }
          profIndex_++;
          std::advance(it_profidx_sorted, 1);
        }
      } else {
        // This will not work for pressures in ascending order
        throw eckit::BadParameter("sort order is ascending.", Here());
      }
    }

    // Number of levels to which QC checks should be applied
    numProfileLevels_ = static_cast<int> (profileIndices_.size());

    if (numProfileLevels_ > 0) {
      oops::Log::debug() << "First and last profile indices: " << profileIndices_.front()
                         << ", " << profileIndices_.back() << std::endl;
    }

    // Replace with maxlev if defined (a legacy of the OPS code)
    if (options_.maxlev.value() != boost::none) {
      numProfileLevels_ = std::min(options_.maxlev.value().get(), numProfileLevels_);
    }

    // Update counters and iterators (if used)

    // Current profile number
    profileNumCurrent_  = profileNumToFind_;

    // Return value indicates whether or not the end of the entire sample has been reached
    if (profIndex_ < profileNums_.size()) {
      // Next profile number to find
      profileNumToFind_ = profileNums_[profIndex_];
      // Iterator corresponding to next profile
      if (!obsdb_.obs_sort_var().empty() &&
          obsdb_.obs_sort_order() == "descending") {
        std::advance(profidx_current_, 1);
      }
    }
  }

  void ProfileIndices::validateTotalNumProf() {
    // If no sorting is performed on the observations it is possible that
    // two separate profiles could (inadvertently) have the same value of the group variable.
    // This would lead to an incorrect value of nrecs.
    // This routine ensures that nrecs is the same as the actual number of profiles in the sample.

    size_t profNum = 0;
    std::vector <size_t> allProfileNums = {0};
    for (size_t j = 0; j < profileNums_.size(); ++j) {
      size_t profNum_j = profileNums_[j];
      if (profNum_j != profNum) {
        allProfileNums.emplace_back(profNum_j);
        profNum = profNum_j;
      }
    }

    if (allProfileNums.size() != obsdb_.nrecs()) {
      throw eckit::BadValue("incorrect number of profiles in unsorted sample", Here());
    }
  }
}  // namespace ufo
