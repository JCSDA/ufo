/*
 * (C) Copyright 2021 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_UTILS_PRIMITIVEVARIABLES_H_
#define UFO_UTILS_PRIMITIVEVARIABLES_H_

#include <memory>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/ArrowProxy.h"

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {

class ObsFilterData;

/// A proxy object providing access to the name, group and values of a "primitive", i.e.
/// single-channel, variable pointed to by a PrimitiveVariablesIterator.
///
/// \see PrimitiveVariables for an example of its use.
///
/// \warning Never use a PrimitiveVariable object obtained by dereferencing a
/// PrimitiveVariablesIterator after incrementing this iterator or when the iterator has gone out
/// of scope.
class PrimitiveVariable {
 public:
  PrimitiveVariable(const Variable &variable, size_t channelIndex,
                    const std::vector<float> &values)
    : variable_(variable), channelIndex_(channelIndex), values_(values)
  {}

  /// Return the name of the primitive variable.
  std::string name() const { return variable_.variable(channelIndex_); }
  /// Return the group of the primitive variable.
  const std::string &group() const { return variable_.group(); }
  /// Return a Variable object representing the primitive variable.
  Variable variable() const { return variable_[channelIndex_]; }
  /// Return the values of the primitive variable at all observation locations.
  const std::vector<float> &values() const { return values_; }

 private:
  const Variable &variable_;
  size_t channelIndex_;
  const std::vector<float> &values_;
};

/// \brief Iterator over the names and values of primitive variables held in a Variables object.
///
/// \see PrimitiveVariables
///
/// Note: this iterator exists to support range-based for loops, but it doesn't have a copy
/// constructor (it would need to be very costly -- requiring a copy of the ObsDataVector stored
/// in the iterator), so it can't be passed to STL algorithms that take iterators by value.
class PrimitiveVariablesIterator {
 public:
  typedef ptrdiff_t difference_type;
  typedef PrimitiveVariable value_type;
  typedef PrimitiveVariable reference;
  typedef ArrowProxy<PrimitiveVariable> pointer;
  typedef std::forward_iterator_tag iterator_category;

  struct BeginTag {};
  struct EndTag {};

  /// \brief Create an iterator pointing to the first primitive variable in \p variables.
  PrimitiveVariablesIterator(const Variables &variables, const ObsFilterData &data, BeginTag)
    : variables_(variables), data_(data), variableIndex_(0), channelIndex_(0)
  {
    loadCurrentVariable();
  }

  /// \brief Create an iterator pointing past the range of primitive variables in \p variables.
  PrimitiveVariablesIterator(const Variables &variables, const ObsFilterData &data, EndTag)
     : variables_(variables), data_(data), variableIndex_(variables_.size()), channelIndex_(0)
  {}

  /// \brief Dereference the iterator, returning a proxy object whose methods such as name(),
  /// group() and values() can be called to get access to the name, group and values of the
  /// primitive variable pointed to by the iterator.
  PrimitiveVariable operator*() const {
    return PrimitiveVariable(variables_[variableIndex_], channelIndex_, (*vector_)[channelIndex_]);
  }

  ArrowProxy<PrimitiveVariable> operator->() const {
    return ArrowProxy<PrimitiveVariable>(operator*());
  }

  PrimitiveVariablesIterator& operator++() {
    if (variableIndex_ < variables_.size()) {
      ++channelIndex_;
      if (channelIndex_ == variables_[variableIndex_].size()) {
        ++variableIndex_;
        channelIndex_ = 0;
        loadCurrentVariable();
      }
    }
    return *this;
  }

  bool operator==(const PrimitiveVariablesIterator& other) const {
    // In principle we should also check variables_ and data_ for equality.
    // We don't do it for efficiency.
    return variableIndex_ == other.variableIndex_;
    return channelIndex_ == other.channelIndex_;
  }

  bool operator!=(const PrimitiveVariablesIterator& other) const {
    return !operator==(other);
  }

 private:
  void loadCurrentVariable();

 private:
  const Variables &variables_;
  const ObsFilterData &data_;
  size_t variableIndex_;
  size_t channelIndex_;
  std::unique_ptr<ioda::ObsDataVector<float>> vector_;
};

/// \brief A range covering all "primitive" (single-channel) variables in a Variables object.
///
/// It makes it possible to iterate over all these primitive variables using a range-based for loop.
/// Each iteration has access to the name and group of a single variable as well as the vector of
/// its values at all observation locations. For example:
///
///     Variables vars = ...;
///     ObsFilterData data = ...;  // within a filter, we'd typically use the data_ member variable
///     for (PrimitiveVariable var : PrimitiveVariables(vars, data)) {
///       // Access to the name of the current primitive variable: var.name()
///       // Access to its group: var.group()
///       // Access to its values (a vector of floats): var.values()
///     }
class PrimitiveVariables {
 public:
  PrimitiveVariables(const Variables &variables, const ObsFilterData &data)
    : variables_(variables), data_(data)
  {}

  PrimitiveVariablesIterator begin() const {
    return PrimitiveVariablesIterator(variables_, data_, PrimitiveVariablesIterator::BeginTag());
  }

  PrimitiveVariablesIterator end() const {
    return PrimitiveVariablesIterator(variables_, data_, PrimitiveVariablesIterator::EndTag());
  }

 private:
  const Variables &variables_;
  const ObsFilterData &data_;
};

}  // namespace ufo

#endif  // UFO_UTILS_PRIMITIVEVARIABLES_H_
