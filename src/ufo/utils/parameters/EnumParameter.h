/*
 * (C) Copyright 2019 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_UTILS_PARAMETERS_ENUMPARAMETER_H_
#define UFO_UTILS_PARAMETERS_ENUMPARAMETER_H_

#include <string>

#include "eckit/config/Configuration.h"
#include "ufo/utils/parameters/ParameterBase.h"

namespace ufo {

/// \brief Convert a string to a member of the enumeration type \c EnumType.
///
/// This template must be specialised for all enumeration types for which EnumParameter is
/// expected to work.
template <typename EnumType>
EnumType enumFromString(const std::string &s);

/// \brief An enumeration-valued option.
template <typename EnumType>
class EnumParameter : public ParameterBase {
 public:
  EnumParameter(const char *name, const EnumType& defaultValue, Parameters *parent = nullptr)
    : ParameterBase(parent), name_(name), value_(defaultValue)
  {}

  void deserialize(const eckit::Configuration &config) override {
    std::string value;
    if (config.get(name_, value)) {
      value_ = enumFromString<EnumType>(value);
    }
  }

  EnumType value() const { return value_; }

  operator EnumType() const { return value_; }

 private:
  std::string name_;
  EnumType value_;
};

}  // namespace ufo

#endif  // UFO_UTILS_PARAMETERS_ENUMPARAMETER_H_
