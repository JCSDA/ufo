/*
 * (C) Copyright 2019 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/utils/RecursiveSplitter.h"

#include <algorithm>
#include <numeric>

namespace ufo
{

RecursiveSplitter::RecursiveSplitter(size_t numIds) {
  orderedIds_.resize(numIds);
  std::iota(orderedIds_.begin(), orderedIds_.end(), 0);
  initializeEncodedGroups();
}

void RecursiveSplitter::initializeEncodedGroups() {
  const size_t numIds = orderedIds_.size();
  encodedGroups_.resize(numIds);
  if (numIds == 1) {
    // There are no groups of 2 or more elements
    encodedGroups_[0] = numIds;
  } else if (numIds > 1) {
    encodedGroups_[0] = 0;
    encodedGroups_[1] = numIds - 1;
  }
}

template <typename T>
void RecursiveSplitter::groupByImpl(const std::vector<T> &categories) {
  auto orderedCategory = [&](size_t index) { return categories[orderedIds_[index]]; };

  const auto numIds = orderedIds_.size();
  ptrdiff_t lastIndexInGroup = -1;
  ptrdiff_t lastIndexInLastGroup = -1;
  while (lastIndexInGroup + 1 < numIds) {
    const size_t firstIndexInGroup = encodedGroups_[lastIndexInGroup + 1];
    if (firstIndexInGroup + 1 >= numIds)
      break;
    lastIndexInGroup = encodedGroups_[firstIndexInGroup + 1];

    std::sort(orderedIds_.begin() + firstIndexInGroup,
              orderedIds_.begin() + lastIndexInGroup + 1,
              [&categories](size_t idA, size_t idB) { return categories[idA] < categories[idB];});

    // Now update the group
    size_t newFirstIndex = firstIndexInGroup;
    for (size_t newLastIndex = firstIndexInGroup;
         newLastIndex <= lastIndexInGroup;
         ++newLastIndex) {
      if (newLastIndex != lastIndexInGroup &&
          orderedCategory(newLastIndex) == orderedCategory(newLastIndex + 1))
          continue;

      if (newLastIndex > newFirstIndex) {
        // We've found a new group of equivalent elements
        encodedGroups_[lastIndexInLastGroup + 1] = newFirstIndex;
        encodedGroups_[newFirstIndex + 1] = newLastIndex;
        lastIndexInLastGroup = newLastIndex;
      }
      newFirstIndex = newLastIndex + 1;
    }
  }

  if (lastIndexInLastGroup + 1 < numIds)
    encodedGroups_[lastIndexInLastGroup + 1] = numIds;
}

template void RecursiveSplitter::groupByImpl(const std::vector<int> &);
template void RecursiveSplitter::groupByImpl(const std::vector<size_t> &);

}  // namespace ufo
