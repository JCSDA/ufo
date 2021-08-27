/*
 * (C) Crown copyright 2021, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_UTILS_METOFFICE_METOFFICESORT_H_
#define UFO_UTILS_METOFFICE_METOFFICESORT_H_

#include <utility>

namespace ufo {

namespace metofficesortdetail {

/// Sink element `parentIndex` of the heap starting at `heapStart` with `heapLength` elements
/// until its correct position.
///
/// If the parent's key is smaller than the key of any of its children, the parent is swapped with
/// its largest child. This is repeated until the parent's key is at least as large as the keys of
/// its children.
template <typename RandomIt, typename UnaryOperation>
void sinkParent(typename RandomIt::difference_type parentIndex,
                RandomIt heapStart, typename RandomIt::difference_type heapLength,
                const UnaryOperation &key) {
  typedef typename RandomIt::difference_type Index;

  const auto parentKey = key(*(heapStart + parentIndex));
  Index childIndex = 2 * parentIndex + 1;  // left child
  while (childIndex < heapLength) {
    // Identify the largest child (preferring the left child if they're equal)
    auto childKey = key(*(heapStart + childIndex));
    Index rightChildIndex = childIndex + 1;
    if (rightChildIndex < heapLength) {
      auto rightChildKey = key(*(heapStart + rightChildIndex));
      if (rightChildKey > childKey) {
        childIndex = rightChildIndex;
        childKey = std::move(rightChildKey);
      }
    }
    // Is the parent at least as large as its largest child?
    if (parentKey >= childKey)
      break;  // Yes. The parent is now at the correct position; we're done
    // No. Swap the parent with its largest child
    std::swap(*(heapStart + parentIndex), *(heapStart + childIndex));
    // Continue sinking the parent.
    parentIndex = childIndex;
    childIndex = 2 * parentIndex + 1;
  }
}

/// Arrange the range `[first, last)` into a max-heap ordered by key \p key.
template <typename RandomIt, typename UnaryOperation>
void makeHeap(RandomIt first, RandomIt last, const UnaryOperation &key) {
  const auto len = last - first;
  // Arrange all trees with more than one level into heaps
  for (auto parentIndex = len / 2 - 1; parentIndex >= 0; --parentIndex) {
    sinkParent(parentIndex, first, len, key);
  }
}

/// Sort the heap `[first, last)` in the order of ascending keys
template <typename RandomIt, typename UnaryOperation>
void sortHeap(RandomIt first, RandomIt last, const UnaryOperation &key) {
  auto len = last - first;
  while (len > 1) {
    // Move the largest element to the end of the heap and exclude it from the heap.
    std::swap(*first, *(first + --len));
    // Arrange the remaining elements into a heap again.
    sinkParent(0, first, len, key);
  }
}

}  // namespace metofficesortdetail

/// \brief Sort the range `[first, last)` in the order of ascending keys using the same algorithm as
/// the `Ops_Integer/Char/RealSort` subroutines from the Met Office OPS system.
///
/// \param first, last
///   The range of elements to sort.
/// \param key
///   An unary functor taking an element of the range (a dereferenced iterator `RandomIt`)
///   and returning a key used as the sorting criterion.
template <typename RandomIt, typename UnaryOperation>
void metOfficeSort(RandomIt first, RandomIt last, const UnaryOperation &key) {
  metofficesortdetail::makeHeap(first, last, key);
  metofficesortdetail::sortHeap(first, last, key);
}

/// \brief Sort the range `[first, last)` in ascending order using the same algorithm as
/// the `Ops_Integer/Char/RealSort` subroutines from the Met Office OPS system.
///
/// \param first, last
///   The range of elements to sort.
template <typename RandomIt>
void metOfficeSort(RandomIt first, RandomIt last) {
  metOfficeSort(first, last, [](const auto &x) { return x; });
}

}  // namespace ufo

#endif  // UFO_UTILS_METOFFICE_METOFFICESORT_H_
