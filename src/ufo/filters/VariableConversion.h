/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_VARIABLECONVERSION_H_
#define UFO_FILTERS_VARIABLECONVERSION_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"

namespace eckit {
class Configuration;
}

namespace ioda {
template <typename DATATYPE>
class ObsDataVector;
class ObsSpace;
}

namespace ufo {
class VariableConversionParameters;
}

namespace ufo {

/// \brief Main filter to apply some variable convertion.
///
/// See VariableConversionParameters for the documentation of the available
/// parameters and options.
///
/// \par Important:
/// Any new variable created is assigned to the observation space with the
/// "@DerivedValue" tag.
///
class VariableConversion : public FilterBase,
                           private util::ObjectCounter<VariableConversion> {
 public:
  static const std::string classname() { return "ufo::VariableConversion"; }
  // This Constructor function initializes an instance of the
  // filter based on options specified in the YAML configuration file.
  VariableConversion(ioda::ObsSpace &, const eckit::Configuration &,
                     std::shared_ptr<ioda::ObsDataVector<int>>,
                     std::shared_ptr<ioda::ObsDataVector<float>>);
  // Destructor
  ~VariableConversion();

 private:
  /// Configurable options
  std::unique_ptr<VariableConversionParameters> options_;
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override { return QCflags::pass; }
};

}  // namespace ufo

#endif  // UFO_FILTERS_VARIABLECONVERSION_H_
