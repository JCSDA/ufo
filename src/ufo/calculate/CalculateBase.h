/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_CALCULATE_CALCULATEBASE_H_
#define UFO_CALCULATE_CALCULATEBASE_H_

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

#include "ufo/calculate/Formulas.h"
#include "ufo/filters/QCflags.h"
#include "ufo/filters/VariableConversionParameters.h"

namespace ioda {
template <typename DATATYPE>
class ObsDataVector;
class ObsSpace;
class ObsVector;
}

namespace ufo {
class VariableConversionParameters;
}

namespace ufo {

/// \brief Base class for variable conversion
class CalculateBase {
 public:
  CalculateBase(const VariableConversionParameters &options,
                ioda::ObsSpace &os,
                const std::shared_ptr<ioda::ObsDataVector<int>>& flags);
  /// Destructor
  virtual ~CalculateBase() {}
  /// Run variable conversion
  virtual void runCalculate() = 0;

 private:
  void filterObservation(const std::string &variable, std::vector<float> &obsVector) const;
  /// Method used for the calculation
  formulas::Method method_;
  /// The observation name
  std::string obsName_;

 protected:  // variables
  void getObservation(const std::string &originalTag, const std::string &varName,
                      std::vector<float> &obsVector, bool require = false) const;
  /// subclasses to access Method used for the calculation
  formulas::Method method() const { return method_; }
  /// subclasses to access the observation name
  std::string obsName() const { return obsName_; }
  /// Configurable parameters
  const VariableConversionParameters &options_;
  /// Observation space
  ioda::ObsSpace &obsdb_;
  const ioda::ObsDataVector<int> &flags_;
  /// Missing value (float)
  const float missingValueFloat = util::missingValue(1.0f);
  /// output tag for derived parameters
  const std::string outputTag = "DerivedValue";
};

/// \brief Calculate factory
class CalculateFactory {
 public:
  static std::unique_ptr<CalculateBase> create(
      const std::string &, const VariableConversionParameters &,
      ioda::ObsSpace &, const std::shared_ptr<ioda::ObsDataVector<int>> &);
  virtual ~CalculateFactory() = default;

 protected:
  explicit CalculateFactory(const std::string &);

 private:
  virtual std::unique_ptr<CalculateBase> make(
      const VariableConversionParameters &, ioda::ObsSpace &,
      const std::shared_ptr<ioda::ObsDataVector<int>> &) = 0;
  static std::map<std::string, CalculateFactory *> &getMakers() {
    static std::map<std::string, CalculateFactory *> makers_;
    return makers_;
  }
};

/// \brief Calculate maker
template <class T>
class CalculateMaker : public CalculateFactory {
  virtual std::unique_ptr<CalculateBase> make(
      const VariableConversionParameters &options, ioda::ObsSpace &os,
      const std::shared_ptr<ioda::ObsDataVector<int>> &flags) {
    return std::unique_ptr<CalculateBase>(new T(options, os, flags));
  }

 public:
  explicit CalculateMaker(const std::string &name) : CalculateFactory(name) {}
};

}  // namespace ufo

#endif  // UFO_CALCULATE_CALCULATEBASE_H_
