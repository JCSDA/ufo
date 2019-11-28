/*
 * (C) Copyright 2019 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_UTILS_PARAMETERS_OPTIONALPARAMETERVARIABLE_H_
#define UFO_UTILS_PARAMETERS_OPTIONALPARAMETERVARIABLE_H_

#include <string>

#include "oops/util/parameters/OptionalParameter.h"
#include "ufo/filters/Variable.h"

namespace oops {

template <>
class OptionalParameter<ufo::Variable> : public ParameterBase {
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

  boost::optional<ufo::Variable> value() const { return value_; }

  operator boost::optional<ufo::Variable>() const { return value_; }

 private:
  std::string name_;
  boost::optional<ufo::Variable> value_;
};

}  // namespace oops

#endif  // UFO_UTILS_PARAMETERS_OPTIONALPARAMETERVARIABLE_H_
