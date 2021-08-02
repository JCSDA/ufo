/*
 * (C) Crown copyright 2021, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_UTILS_METOFFICE_METOFFICESORT_H_
#define UFO_UTILS_METOFFICE_METOFFICESORT_H_

#include <algorithm>

namespace ufo {

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
  auto compare = [&key] (const auto &a, const auto &b) { return key(a) <= key(b); };
  std::make_heap(first, last, compare);
  std::sort_heap(first, last, compare);
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
