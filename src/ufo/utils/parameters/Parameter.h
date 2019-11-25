/*
 * (C) Copyright 2019 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_UTILS_PARAMETERS_PARAMETER_H_
#define UFO_UTILS_PARAMETERS_PARAMETER_H_

#include <string>

#include "eckit/config/Configuration.h"
#include "ufo/utils/parameters/ParameterBase.h"

namespace ufo {

/// \brief A parameter with a default value.
template <typename T>
class Parameter : public ParameterBase {
 public:
  Parameter(const char *name, const T& defaultValue, Parameters *parent = nullptr)
    : ParameterBase(parent), name_(name), value_(defaultValue)
  {}

  void deserialize(const eckit::Configuration &config) override {
    config.get(name_, value_);
  }

  T value() const { return value_; }

  operator T() const { return value_; }

 private:
  std::string name_;
  T value_;
};

}  // namespace ufo

#endif  // UFO_UTILS_PARAMETERS_PARAMETER_H_
