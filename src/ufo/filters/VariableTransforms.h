/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_VARIABLETRANSFORMS_H_
#define UFO_FILTERS_VARIABLETRANSFORMS_H_

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
class VariableTransformParametersBase;
}

namespace ufo {

/// \brief Main filter to apply some variable conversion.
///
/// See VariableTransformParametersBase for the documentation of the available
/// parameters and options.
///
/// \par Important:
/// Any new variable created is assigned to the observation space with the
/// "DerivedObsValue/" tag.
///
class VariableTransforms : public FilterBase,
                           private util::ObjectCounter<VariableTransforms> {
 public:
  static const std::string classname() { return "ufo::VariableTransforms"; }
  // This Constructor function initializes an instance of the
  // filter based on options specified in the YAML configuration file.
  VariableTransforms(ioda::ObsSpace &, const eckit::Configuration &,
                     std::shared_ptr<ioda::ObsDataVector<int>>,
                     std::shared_ptr<ioda::ObsDataVector<float>>);
  // Destructor
  ~VariableTransforms();

 private:
  std::unique_ptr<VariableTransformParametersBase> parameters_;
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override { return QCflags::pass; }
};

}  // namespace ufo

#endif  // UFO_FILTERS_VARIABLETRANSFORMS_H_
