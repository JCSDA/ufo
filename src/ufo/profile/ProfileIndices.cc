/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <utility>

#include "ufo/profile/ProfileIndices.h"

#include "ufo/utils/Flags.h"

namespace ufo {
  ProfileIndices::ProfileIndices(ioda::ObsSpace &obsdb,
                                 const ProfileConsistencyCheckParameters &options,
                                 const std::vector <bool> &apply)
    : obsdb_(obsdb),
      options_(options),
      apply_(apply),
      profileNums_(obsdb.recnum()),
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

    // If sorting observations, ensure indices remain ascending
    if (!obsdb_.obs_sort_var().empty() &&
        obsdb_.obs_sort_order() == "descending") {
      descendingSortWithAscendingIndices();
    }
  }

  void ProfileIndices::determineProfileIndices()
  {
    profileIndices_.clear();

    // If there are no profiles in the sample, warn and exit
    if (profileNums_.size() == 0) {
      oops::Log::warning() << "No profiles found in the sample" << std::endl;
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
        throw eckit::NotImplemented("ascending profile sort order", Here());
      }
    }

    // Number of levels to which QC checks should be applied
    numLevelsToCheck_ = static_cast<int> (profileIndices_.size());

    // Replace with maxlev if defined (a legacy of the OPS code)
    if (options_.maxlev.value() != boost::none) {
      numLevelsToCheck_ = std::min(options_.maxlev.value().get(), numLevelsToCheck_);
    }

    // Update counters and iterators (if used)
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

  void ProfileIndices::descendingSortWithAscendingIndices() {
    // This is a modified version of ioda::ObsData::BuildSortedObsGroups.
    // The sort variable can in theory be anything but it is usually air_pressure.
    // The code sorts such that pressure is descending
    // but, when at least two pressures are equal,
    // the associated profile indices are sorted to be ascending
    // (unlike the version in BuildSortedObsGroups for which the indices descend).

    const size_t nlocs = obsdb_.nlocs();

    // Get the sort variable from the data store, and convert to a vector of floats.
    std::vector<float> SortValues(nlocs);
    obsdb_.get_db("MetaData", obsdb_.obs_sort_var(), SortValues);

    // Construct a temporary structure to do the sorting, then transfer the results
    // to the data member profidx_.
    TmpProfIdxMap TmpProfIdx;
    for (size_t iloc = 0; iloc < nlocs; iloc++) {
      TmpProfIdx[profileNums_[iloc]].push_back(std::make_pair(SortValues[iloc], iloc));
    }

    for (TmpProfIdxIter iprof = TmpProfIdx.begin(); iprof != TmpProfIdx.end(); ++iprof) {
      // Use a lambda function to access implement a descending order sort
      // In BuildSortedObsGroups this is just p1 > p2
      // This implementation reverses the treatment of the indices
      // i.e. p1.second < p2.second would usually be p1.second > p2.second
      sort(iprof->second.begin(), iprof->second.end(),
           [](const std::pair<float, std::size_t> & p1,
              const std::pair<float, std::size_t> & p2) {
             return(p2.first < p1.first ||
                    (!(p2.first < p1.first) && p1.second < p2.second));});
    }

    // Copy indexing to the profidx_ data member.
    for (TmpProfIdxIter iprof = TmpProfIdx.begin(); iprof != TmpProfIdx.end(); ++iprof) {
      profidx_[iprof->first].resize(iprof->second.size());
      for (std::size_t iloc = 0; iloc < iprof->second.size(); iloc++) {
        profidx_[iprof->first][iloc] = iprof->second[iloc].second;
      }
    }

    // Iterator pointing to beginning of indices
    profidx_current_ = profidx_.begin();
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
