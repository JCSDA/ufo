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
                const std::shared_ptr<ioda::ObsDataVector<int>>& flags,
                const std::vector<bool>& apply);
  /// Destructor
  virtual ~TransformBase() {}
  /// Run variable conversion
  virtual void runTransform() = 0;

 private:
  /// Method used for the calculation
  formulas::MethodFormulation method_;
  formulas::MethodFormulation formulation_;
  bool UseValidDataOnly_;
  bool AllowSuperSaturation_;
  /// The observation name
  std::string obsName_;
  /// templated function for float, int data types
  template <typename T>
  void filterObservation(const std::string &variableName,
                         std::vector<T> &obsVector) const {
    if (flags_.has(variableName)) {
      constexpr T one = T(1.0);
      const T missing = util::missingValue(one);
      const std::vector<int> *varFlags = &flags_[variableName];

      std::transform(obsVector.begin(), obsVector.end(),  // Input range 1
                     varFlags->begin(),  // 1st element of input range vector 2 (must be same size)
                     obsVector.begin(),  // 1st element of output range (must be same size)
                     [missing](T obsvalue, int flag)
                     { return flag == QCflags::missing || flag == QCflags::bounds
                       ? missing : obsvalue; });
    }
  }

 protected:  // variables
  /// subclasses to access Method and formualtion used for the calculation
  formulas::MethodFormulation method() const { return method_; }
  formulas::MethodFormulation formulation() const { return formulation_; }
  bool UseValidDataOnly() const { return UseValidDataOnly_; }
  bool AllowSuperSaturation() const { return AllowSuperSaturation_; }
  void SetUseValidDataOnly(bool t) {UseValidDataOnly_ = t; }
  /// subclasses to access the observation name
  std::string obsName() const { return obsName_; }
  /// Configurable parameters
  const VariableTransformsParameters &options_;
  /// Observation space
  ioda::ObsSpace &obsdb_;
  const ioda::ObsDataVector<int> &flags_;
  const std::vector<bool> &apply_;
  /// Missing value (int)
  const int missingValueInt = util::missingValue(1);
  /// Missing value (float)
  const float missingValueFloat = util::missingValue(1.0f);
  /// output tag for derived parameters
  const std::string outputTag = "DerivedValue";
  /// templated function for float, int data types
  template <typename T>
  void getObservation(const std::string &originalTag, const std::string &varName,
                      std::vector<T> &obsVector, bool require = false) const {
    const size_t nlocs = obsdb_.nlocs();

    if (obsdb_.has("DerivedValue", varName)) {
      obsVector = std::vector<T>(nlocs);
      obsdb_.get_db("DerivedValue", varName, obsVector);
      // Set obsValue to missingValueFloat if flag is equal to QCflags::missing or QCflags::bounds
      if (UseValidDataOnly()) filterObservation(varName, obsVector);
    } else if (obsdb_.has(originalTag, varName)) {
      obsVector = std::vector<T>(nlocs);
      obsdb_.get_db(originalTag, varName, obsVector);
      // Set obsValue to missingValueFloat if flag is equal to QCflags::missing or QCflags::bounds
      if (UseValidDataOnly()) filterObservation(varName, obsVector);
    }

    if (require && obsVector.empty()) {
      throw eckit::BadValue("The parameter `" + varName + "@" + originalTag +
                          "` does not exist in the ObsSpace ", Here());
    }
  }
};

/// \brief Transform factory
class TransformFactory {
 public:
  static std::unique_ptr<TransformBase> create(
      const std::string &, const VariableTransformsParameters &,
      ioda::ObsSpace &, const std::shared_ptr<ioda::ObsDataVector<int>> &,
      const std::vector<bool> &);
  virtual ~TransformFactory() = default;

 protected:
  explicit TransformFactory(const std::string &);

 private:
  virtual std::unique_ptr<TransformBase> make(
      const VariableTransformsParameters &, ioda::ObsSpace &,
      const std::shared_ptr<ioda::ObsDataVector<int>> &,
      const std::vector<bool> &) = 0;
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
      const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
      const std::vector<bool> &apply) {
    return std::unique_ptr<TransformBase>(new T(options, os, flags, apply));
  }

 public:
  explicit TransformMaker(const std::string &name) : TransformFactory(name) {}
};

}  // namespace ufo

#endif  // UFO_VARIABLETRANSFORMS_TRANSFORMBASE_H_
