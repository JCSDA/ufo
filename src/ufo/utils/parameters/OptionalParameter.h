/*
 * (C) Copyright 2019 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_UTILS_PARAMETERS_OPTIONALPARAMETER_H_
#define UFO_UTILS_PARAMETERS_OPTIONALPARAMETER_H_


#include <string>

#include <boost/optional.hpp>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/parameters/ParameterBase.h"

namespace ufo {

/// \brief A parameter that does not need to be present in a Configuration and for which no
/// sensible default value can be defined.
template <typename T>
class OptionalParameter : public ParameterBase {
 public:
  explicit OptionalParameter(const char *name, Parameters *parent = nullptr)
    : ParameterBase(parent), name_(name)
  {}

  void deserialize(const eckit::Configuration &config) override {
    T value;
    if (config.get(name_, value)) {
      value_ = value;
    }
  }

  boost::optional<T> value() const { return value_; }

  operator boost::optional<T>() const { return value_; }

 private:
  std::string name_;
  boost::optional<T> value_;
};

template <>
class OptionalParameter<util::Duration> : public ParameterBase {
 public:
  explicit OptionalParameter(const char *name, Parameters *parent = nullptr)
    : ParameterBase(parent), name_(name)
  {}

  void deserialize(const eckit::Configuration &config) override {
    std::string value;
    if (config.get(name_, value)) {
      value_ = util::Duration(value);
    }
  }

  boost::optional<util::Duration> value() const { return value_; }

  operator boost::optional<util::Duration>() const { return value_; }

 private:
  std::string name_;
  boost::optional<util::Duration> value_;
};

template <>
class OptionalParameter<util::DateTime> : public ParameterBase {
 public:
  explicit OptionalParameter(const char *name, Parameters *parent = nullptr)
    : ParameterBase(parent), name_(name)
  {}

  void deserialize(const eckit::Configuration &config) override {
    std::string value;
    if (config.get(name_, value)) {
      value_ = util::DateTime(value);
    }
  }

  boost::optional<util::DateTime> value() const { return value_; }

  operator boost::optional<util::DateTime>() const { return value_; }

 private:
  std::string name_;
  boost::optional<util::DateTime> value_;
};

template <>
class OptionalParameter<Variable> : public ParameterBase {
 public:
  explicit OptionalParameter(const char *name, Parameters *parent = nullptr)
    : ParameterBase(parent), name_(name)
  {}

  void deserialize(const eckit::Configuration &config) override {
    if (config.has(name_)) {
      eckit::LocalConfiguration varConf(config, name_);
      if (!varConf.has("name")) {
        // TODO(wsmigaj): shouldn't ufo::Variable itself throw an exception if
        // the 'name' property is not specified?
        throw eckit::BadParameter("No variable name specified", Here());
      }
      ufo::Variable var(varConf);
      value_ = var;
    }
  }

  boost::optional<Variable> value() const { return value_; }

  operator boost::optional<Variable>() const { return value_; }

 private:
  std::string name_;
  boost::optional<Variable> value_;
};

}  // namespace ufo

#endif  // UFO_UTILS_PARAMETERS_OPTIONALPARAMETER_H_
