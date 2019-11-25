/*
 * (C) Copyright 2019 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_UTILS_PARAMETERS_PARAMETERS_H_
#define UFO_UTILS_PARAMETERS_PARAMETERS_H_

#include <vector>

#include "ufo/utils/parameters/ParameterBase.h"

namespace eckit {
  class Configuration;
}

namespace ufo {

/// \brief A collection of parameters that can be loaded together from a Configuration object.
class Parameters : public ParameterBase {
 public:
  /// \brief Add \p parameter to the list of parameters processed by subsequent calls to
  /// deserialize().
  void registerChild(ParameterBase &parameter);

  /// \brief Load all previously registered parameters from \p config.
  void deserialize(const eckit::Configuration &config) override;

 private:
  std::vector<ParameterBase*> children_;
};

}  // namespace ufo

#endif  // UFO_UTILS_PARAMETERS_PARAMETERS_H_
