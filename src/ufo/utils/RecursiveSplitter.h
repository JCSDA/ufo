/*
 * (C) Copyright 2019 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_RECURSIVESPLITTER_H_
#define UFO_UTILS_RECURSIVESPLITTER_H_

#include <algorithm>
#include <cassert>
#include <cstddef>  // for size_t
#include <string>
#include <vector>

#include "ufo/utils/ArrowProxy.h"

namespace ufo {

/// \brief Partitions an array into groups of elements equivalent according to certain criteria.
///
/// Example:
/// \code
/// RecursiveSplitter splitter(8);
/// // Tell the splitter that according to a certain criterion, some elements belong to category #1
/// // and others to category #2
/// splitter.groupBy(std::vector<int>({1, 2, 1, 2, 1, 2, 1, 2}));
/// std::cout << "After first split:" << std::endl;
/// for (const auto &group : splitter.groups()) {
///   std::cout << "  Elements with these indices are equivalent: ";
///   for (const auto &index : group)
///     std::cout << index << ", ";
///   std::cout << std::endl;
/// }
/// // Tell the splitter that according to a different criterion, there is a different split between
/// // categories
/// splitter.groupBy(std::vector<int>({1, 1, 1, 1, 2, 2, 2, 2}));
/// std::cout << "After second split:" << std::endl;
/// for (const auto &group : splitter.groups()) {
///   std::cout << "  Elements with these indices are equivalent: ";
///   for (const auto &index : group)
///     std::cout << index << ", ";
///   std::cout << std::endl;
/// }
/// \endcode
///
/// Output:
/// \code{.unparsed}
/// After first split:
///   Elements with these indices are equivalent: 0, 2, 4, 6,
///   Elements with these indices are equivalent: 1, 3, 5, 7,
/// After second split:
///   Elements with these indices are equivalent: 0, 2,
///   Elements with these indices are equivalent: 4, 6,
///   Elements with these indices are equivalent: 1, 3,
///   Elements with these indices are equivalent: 5, 7,
/// \endcode
///
/// \internal In the implementation, indices into the partitioned array are referred to as _ids_.
/// The term _index_ denotes an index into the vector \c orderedIds_.
class RecursiveSplitter
{
 public:
  /// \brief A range of indices of all array elements belonging to a particular equivalence class.
  class Group
  {
   public:
    Group(const std::vector<size_t>::const_iterator &begin,
          const std::vector<size_t>::const_iterator &end)
      : begin_(begin), end_(end)
    {}

    std::vector<size_t>::const_iterator begin() const
    {
      return begin_;
    }

    std::vector<size_t>::const_iterator end() const
    {
      return end_;
    }

   private:
    std::vector<size_t>::const_iterator begin_;
    std::vector<size_t>::const_iterator end_;
  };

  /// \brief An iterator over equivalence classes consisting of more than one element.
  class MultiElementGroupIterator
  {
   public:
    typedef ptrdiff_t difference_type;
    typedef Group value_type;
    typedef ArrowProxy<Group> reference;
    typedef std::forward_iterator_tag iterator_category;

    struct BeginTag {};
    struct EndTag {};

    explicit MultiElementGroupIterator(const RecursiveSplitter &splitter, BeginTag)
      : splitter_(splitter) {
      if (splitter_.encodedGroups_.empty()) {
        firstIndexInGroup_ = 0;
      } else {
        firstIndexInGroup_ = splitter_.encodedGroups_[0];
      }
    }

    explicit MultiElementGroupIterator(const RecursiveSplitter &splitter, EndTag)
      : splitter_(splitter) {
      firstIndexInGroup_ = splitter_.encodedGroups_.size();
    }

    Group operator*() const {
      assert(!isSentinel());
      size_t lastIndexInGroup = splitter_.encodedGroups_[firstIndexInGroup_ + 1];
      return Group(splitter_.orderedIds_.begin() + firstIndexInGroup_,
                   splitter_.orderedIds_.begin() + lastIndexInGroup + 1);
    }

    ArrowProxy<Group> operator->() const {
      return ArrowProxy<Group>(operator*());
    }

    MultiElementGroupIterator& operator++() {
      const size_t lastIndexInGroup = splitter_.encodedGroups_[firstIndexInGroup_ + 1];
      if (lastIndexInGroup + 1 < splitter_.encodedGroups_.size()) {
        firstIndexInGroup_ = splitter_.encodedGroups_[lastIndexInGroup + 1];
      } else {
        firstIndexInGroup_ = splitter_.encodedGroups_.size();
      }
      return *this;
    }

    bool operator==(const MultiElementGroupIterator& other) const {
      // In principle we should also check splitters for equality. We don't do it for efficiency.
      return firstIndexInGroup_ == other.firstIndexInGroup_;
    }

    bool operator!=(const MultiElementGroupIterator& other) const {
      return !operator==(other);
    }

   protected:
    bool isSentinel() const {
      return firstIndexInGroup_ == splitter_.encodedGroups_.size();
    }

   protected:
    const RecursiveSplitter& splitter_;
    size_t firstIndexInGroup_;
  };

  /// \brief An iterator over all equivalence classes.
  class GroupIterator : private MultiElementGroupIterator {
   public:
    typedef MultiElementGroupIterator::difference_type difference_type;
    typedef MultiElementGroupIterator::value_type value_type;
    typedef MultiElementGroupIterator::reference reference;
    typedef MultiElementGroupIterator::iterator_category iterator_category;

    typedef MultiElementGroupIterator::BeginTag BeginTag;
    typedef MultiElementGroupIterator::EndTag EndTag;

    explicit GroupIterator(const RecursiveSplitter &splitter, BeginTag)
      : MultiElementGroupIterator(splitter, BeginTag()), currentIndex_(0)
    {}

    explicit GroupIterator(const RecursiveSplitter &splitter, EndTag)
      : MultiElementGroupIterator(splitter, EndTag()), currentIndex_(firstIndexInGroup_)
    {}

    Group operator*() const {
      assert(!isSentinel());
      if (currentIndex_ == firstIndexInGroup_) {
        return MultiElementGroupIterator::operator*();
      } else {
        return Group(splitter_.orderedIds_.begin() + currentIndex_,
                     splitter_.orderedIds_.begin() + currentIndex_ + 1);
      }
    }

    ArrowProxy<Group> operator->() const {
      return ArrowProxy<Group>(operator*());
    }

    GroupIterator& operator++() {
      if (currentIndex_ == firstIndexInGroup_) {
        const size_t lastIndexInGroup = splitter_.encodedGroups_[firstIndexInGroup_ + 1];
        if (lastIndexInGroup + 1 < splitter_.encodedGroups_.size()) {
          firstIndexInGroup_ = splitter_.encodedGroups_[lastIndexInGroup + 1];
        } else {
          firstIndexInGroup_ = splitter_.encodedGroups_.size();
        }
        currentIndex_ = lastIndexInGroup + 1;
      } else {
        ++currentIndex_;
      }
      return *this;
    }

    bool operator==(const GroupIterator& other) const {
      // In principle we should also check splitters for equality. We don't do it for efficiency.
      return firstIndexInGroup_ == other.firstIndexInGroup_ &&
             currentIndex_ == other.currentIndex_;
    }

    bool operator!=(const GroupIterator& other) const {
      return !operator==(other);
    }

   private:
    bool isSentinel() const {
      return currentIndex_ == splitter_.encodedGroups_.size();
    }

   private:
    size_t currentIndex_;
  };

  /// \brief A range of equivalence classes.
  class GroupRange
  {
   public:
    explicit GroupRange(const RecursiveSplitter& splitter)
      : splitter_(splitter)
    {}

    GroupIterator begin() const {
      return GroupIterator(splitter_, GroupIterator::BeginTag());
    }

    GroupIterator end() const {
      return GroupIterator(splitter_, GroupIterator::EndTag());
    }
   private:
    const RecursiveSplitter &splitter_;
  };

  /// \brief A range of equivalence classes consisting of more than one element.
  class MultiElementGroupRange
  {
   public:
    explicit MultiElementGroupRange(const RecursiveSplitter& splitter)
      : splitter_(splitter)
    {}

    MultiElementGroupIterator begin() const {
      return MultiElementGroupIterator(splitter_, MultiElementGroupIterator::BeginTag());
    }

    MultiElementGroupIterator end() const {
      return MultiElementGroupIterator(splitter_, MultiElementGroupIterator::EndTag());
    }
   private:
    const RecursiveSplitter &splitter_;
  };

  /// \brief Initialize partitioning of an array of \p numIds elements.
  ///
  /// Initially, all elements are assumed to belong to the same equivalence class.
  explicit RecursiveSplitter(size_t numIds);

  /// \brief Split existing equivalence classes according to a new criterion.
  ///
  /// \param categories
  ///   A vector assigning each element of the partitioned array to a unique category.
  ///
  /// Each existing equivalence class E consisting of N_E elements with indices E_i (0 <= i < N_E)
  /// is split into one or more classes according to the equivalence relation
  ///
  ///   E_i'th element is equivalent to E_j'th element if and only if
  ///   categories[E_i] == categories[E_j].
  void groupBy(const std::vector<size_t> &categories) {
      return groupByImpl(categories);
  }

  /// \overload
  void groupBy(const std::vector<int> &categories) {
      return groupByImpl(categories);
  }

  /// \overload
  void groupBy(const std::vector<std::string> &categories) {
      return groupByImpl(categories);
  }

  /// \brief Return the range of equivalence classes.
  GroupRange groups() const { return GroupRange(*this); }

  /// \brief Return the range of equivalence classes consisting of more than one element.
  MultiElementGroupRange multiElementGroups() const { return MultiElementGroupRange(*this); }

  /// \brief Sort the elements in each equivalence class in ascending order.
  ///
  /// The elements are compared using the binary comparison function \p comp. This function needs
  /// to satisfy the same requirements as the \c comp argument of std::sort().
  template <typename Compare>
  void sortGroupsBy(Compare comp);

  /// \brief Randomly shuffle the elements of each equivalence class.
  void shuffleGroups();

  /// \brief Initialise the random number generator used by shuffleGroups() with a seed.
  ///
  /// \param seed
  ///   Seed with which to initialise the generator.
  /// \param force
  ///   If false, the seed will only be reset if the program has not made any calls to
  ///   util::shuffle() yet.
  void setSeed(unsigned int seed, bool force);

 private:
  void initializeEncodedGroups();

 private:
  template <typename T>
  void groupByImpl(const std::vector<T> &categories);

  /// Indices of elements of the partitioned array ordered by equivalence class.
  std::vector<size_t> orderedIds_;
  /// Encoded locations of multi-element equivalence classes in orderedIds_.
  std::vector<size_t> encodedGroups_;
};

template <typename Compare>
void RecursiveSplitter::sortGroupsBy(Compare comp) {
  for (Group group : multiElementGroups()) {
    std::vector<size_t>::iterator nonConstGroupBegin =
        orderedIds_.begin() + (group.begin() - orderedIds_.cbegin());
    std::vector<size_t>::iterator nonConstGroupEnd =
        orderedIds_.begin() + (group.end() - orderedIds_.cbegin());
    std::sort(nonConstGroupBegin, nonConstGroupEnd, comp);
  }
}

}  // namespace ufo

#endif  // UFO_UTILS_RECURSIVESPLITTER_H_
