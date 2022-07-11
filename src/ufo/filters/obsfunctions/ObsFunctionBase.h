/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONBASE_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONBASE_H_

#include <map>
#include <string>

#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"

namespace util {
  class DateTime;
}

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
}

namespace ufo {
  class Variable;
  class Variables;
  class ObsFilterData;

// -----------------------------------------------------------------------------
/// Base class for computing functions on observation data
///
/// \tparam FunctionValue
///   Type of the values produced by the function. Must be `float`, `int`, `std::string`
///   or `util::DateTime`.
template <typename FunctionValue>
class ObsFunctionBase : private boost::noncopyable {
 public:
  /// Type of the values produced by the function.
  typedef FunctionValue Value_;

  explicit ObsFunctionBase(const eckit::LocalConfiguration conf = eckit::LocalConfiguration()) {}
  virtual ~ObsFunctionBase() {}

/// compute the result of the function
  virtual void compute(const ObsFilterData &,
                       ioda::ObsDataVector<FunctionValue> &) const = 0;

/// geovals required to compute the function
  virtual const ufo::Variables & requiredVariables() const = 0;
};

// -----------------------------------------------------------------------------

/// \brief Common properties of ObsFunctions producing values of type `FunctionValue`.
template <typename FunctionValue>
struct ObsFunctionTraits;

template <>
struct ObsFunctionTraits<float>{
  /// Name of the type of values produced by subclasses of ObsFunctionBase<float>.
  static const char *valueTypeName;
  /// Name of the group identifying ObsFunctions producing floats.
  static const char *groupName;
};

template <>
struct ObsFunctionTraits<int>{
  /// Name of the type of values produced by subclasses of ObsFunctionBase<int>.
  static const char *valueTypeName;
  /// Name of the group identifying ObsFunctions producing ints.
  static const char *groupName;
};

template <>
struct ObsFunctionTraits<std::string>{
  /// Name of the type of values produced by subclasses of ObsFunctionBase<std::string>.
  static const char *valueTypeName;
  /// Name of the group identifying ObsFunctions producing strings.
  static const char *groupName;
};

template <>
struct ObsFunctionTraits<util::DateTime>{
  /// Name of the type of values produced by subclasses of ObsFunctionBase<util::DateTime>.
  static const char *valueTypeName;
  /// Name of the group identifying ObsFunctions producing datetimes.
  static const char *groupName;
};

// -----------------------------------------------------------------------------

/// Factory of ObsFunctions producing values of type `FunctionValue`.
template <typename FunctionValue>
class ObsFunctionFactory {
 public:
  static ObsFunctionBase<FunctionValue> * create(const Variable &);
  virtual ~ObsFunctionFactory() = default;
  static bool functionExists(const std::string &);
 protected:
  explicit ObsFunctionFactory(const std::string &);
 private:
  virtual ObsFunctionBase<FunctionValue> * make(const eckit::LocalConfiguration conf) = 0;
  static std::map < std::string, ObsFunctionFactory * > & getMakers() {
    static std::map < std::string, ObsFunctionFactory * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template <class T>
class ObsFunctionMaker : public ObsFunctionFactory<typename T::Value_> {
  typedef typename T::Value_ Value_;
  typedef ObsFunctionFactory<Value_> Factory_;

  virtual ObsFunctionBase<Value_> * make(const eckit::LocalConfiguration conf)
    { return new T(conf); }
 public:
  explicit ObsFunctionMaker(const std::string & name)
    : Factory_(name) {}
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONBASE_H_
