/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_UTILS_DATAEXTRACTOR_CONSTRAINEDRANGE_H_
#define UFO_UTILS_DATAEXTRACTOR_CONSTRAINEDRANGE_H_

#include <cassert>

namespace ufo {

/// \brief A range of indices.
class ConstrainedRange {
 public:
  /// \brief Create an unconstrained range of size `size`.
  explicit ConstrainedRange(int size = 0) : size_(size) {
    reset();
  }

  /// \brief Return the index of the first element in the range.
  int begin() const { return begin_; }
  /// \brief Return the index of the element past the end of the range.
  int end() const { return end_; }

  /// \brief Return the range length.
  int size() const { return end_ - begin_; }
  /// \brief Return true if the range is empty, false otherwise.
  bool empty() const { return end_ == begin_; }

  /// \brief Constrain the range.
  void constrain(int newBegin, int newEnd) {
    assert(newBegin >= begin_);
    assert(newEnd <= end_);
    begin_ = newBegin;
    end_ = newEnd;
  }

  /// \brief Remove any constraints, resetting the range to its original size.
  void reset() {
    begin_ = 0;
    end_ = size_;
  }

 private:
  int size_;
  int begin_;
  int end_;
};

}  // namespace ufo

#endif  // UFO_UTILS_DATAEXTRACTOR_CONSTRAINEDRANGE_H_
