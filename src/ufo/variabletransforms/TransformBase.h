/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_VARIABLETRANSFORMS_TRANSFORMBASE_H_
#define UFO_VARIABLETRANSFORMS_TRANSFORMBASE_H_

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"
#include "eckit/types/FloatCompare.h"

#include "ioda/ObsDataVector.h"

#include "oops/util/CompareNVectors.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/PropertiesOfNVectors.h"

#include "ufo/filters/QCflags.h"
#include "ufo/filters/VariableTransformsParameters.h"
#include "ufo/variabletransforms/Formulas.h"

namespace ioda {
template <typename DATATYPE>
class ObsDataVector;
class ObsSpace;
class ObsVector;
}

namespace ufo {
class VariableTransformsParameters;
}

namespace ufo {

/// \brief Base class for variable conversion
class TransformBase {
 public:
  TransformBase(const VariableTransformsParameters &options,
                ioda::ObsSpace &os,
                const std::shared_ptr<ioda::ObsDataVector<int>>& flags);
  /// Destructor
  virtual ~TransformBase() {}
  /// Run variable conversion
  virtual void runTransform() = 0;

 private:
  void filterObservation(const std::string &variable, std::vector<float> &obsVector) const;
  /// Method used for the calculation
  formulas::MethodFormulation method_;
  formulas::MethodFormulation formulation_;
  bool UseValidDataOnly_;
  /// The observation name
  std::string obsName_;

 protected:  // variables
  void getObservation(const std::string &originalTag, const std::string &varName,
                      std::vector<float> &obsVector, bool require = false) const;
  /// subclasses to access Method and formualtion used for the calculation
  formulas::MethodFormulation method() const { return method_; }
  formulas::MethodFormulation formulation() const { return formulation_; }
  bool UseValidDataOnly() const { return UseValidDataOnly_; }
  void SetUseValidDataOnly(bool t) {UseValidDataOnly_ = t; }
  /// subclasses to access the observation name
  std::string obsName() const { return obsName_; }
  /// Configurable parameters
  const VariableTransformsParameters &options_;
  /// Observation space
  ioda::ObsSpace &obsdb_;
  const ioda::ObsDataVector<int> &flags_;
  /// Missing value (float)
  const float missingValueFloat = util::missingValue(1.0f);
  /// output tag for derived parameters
  const std::string outputTag = "DerivedValue";
};

/// \brief Transform factory
class TransformFactory {
 public:
  static std::unique_ptr<TransformBase> create(
      const std::string &, const VariableTransformsParameters &,
      ioda::ObsSpace &, const std::shared_ptr<ioda::ObsDataVector<int>> &);
  virtual ~TransformFactory() = default;

 protected:
  explicit TransformFactory(const std::string &);

 private:
  virtual std::unique_ptr<TransformBase> make(
      const VariableTransformsParameters &, ioda::ObsSpace &,
      const std::shared_ptr<ioda::ObsDataVector<int>> &) = 0;
  static std::map<std::string, TransformFactory *> &getMakers() {
    static std::map<std::string, TransformFactory *> makers_;
    return makers_;
  }
};

/// \brief Transform maker
template <class T>
class TransformMaker : public TransformFactory {
  virtual std::unique_ptr<TransformBase> make(
      const VariableTransformsParameters &options, ioda::ObsSpace &os,
      const std::shared_ptr<ioda::ObsDataVector<int>> &flags) {
    return std::unique_ptr<TransformBase>(new T(options, os, flags));
  }

 public:
  explicit TransformMaker(const std::string &name) : TransformFactory(name) {}
};

}  // namespace ufo

#endif  // UFO_VARIABLETRANSFORMS_TRANSFORMBASE_H_
