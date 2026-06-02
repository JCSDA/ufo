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
#include "oops/util/parameters/HasParameters_.h"
#include "oops/util/PropertiesOfNVectors.h"

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/QCflags.h"
#include "ufo/filters/Variables.h"
#include "ufo/filters/VariableTransformParametersBase.h"
#include "ufo/variabletransforms/Formulas.h"


namespace ioda {
template <typename DATATYPE>
class ObsDataVector;
class ObsSpace;
class ObsVector;
}

namespace ufo {
class VariableTransformParametersBase;
class Variables;
}

namespace ufo {

/// \brief Base class for variable conversion
class TransformBase {
 public:
  TransformBase(const VariableTransformParametersBase &options,
                const ObsFilterData& data,
                const std::shared_ptr<ioda::ObsDataVector<int>>& flags,
                const std::shared_ptr<ioda::ObsDataVector<float>>& obserr);
  /// Destructor
  virtual ~TransformBase() {}
  /// Run variable conversion
  virtual void runTransform(const std::vector<bool> &apply) = 0;
  /// Return list of required geovals
  virtual Variables requiredVariables() const { return Variables(); }

 private:
  /// templated function for float, int data types
  template <typename T>
  void filterObservation(const std::string &variableName,
                         std::vector<T> &obsVector) const {
    if (flags_.has(variableName)) {
      const T missing = util::missingValue<T>();
      const std::vector<int> *varFlags = &flags_[variableName];

      std::transform(obsVector.begin(), obsVector.end(),  // Input range 1
                     varFlags->begin(),  // 1st element of input range vector 2 (must be same size)
                     obsVector.begin(),  // 1st element of output range (must be same size)
                     [missing](T obsvalue, int flag)
                     { return flag == QCflags::missing || flag == QCflags::bounds
                       ? missing : obsvalue; });
    }
  }

  /// Method used for the calculation
  formulas::Method method_;
  bool UseValidDataOnly_;
  /// The observation name
  std::string obsName_;

 protected:
  /// templated function for float, int data types
  template <typename T>
  void getObservation(const std::string &originalTag, const std::string &varName,
                      std::vector<T> &obsVector, bool require = false) const {
    if (!obsdb_.has(originalTag, varName)) {
      if (require)
        throw eckit::BadValue("The parameter `" + originalTag + "/" + varName +
                              "` does not exist in the ObsSpace ", Here());
      else
        return;
    }

    obsVector.resize(obsdb_.nlocs());
    obsdb_.get_db(originalTag, varName, obsVector);
    // Set obsValue to missingValue if flag is equal to QCflags::missing or QCflags::bounds
    if (UseValidDataOnly()) filterObservation(varName, obsVector);
  }

  /// \brief Save a transformed variable to the `DerivedObsValue` group of the obs space.
  ///
  /// If the saved variable is a simulated variable, QC flags previously set to `missing` are reset
  /// to `pass` at locations where a valid obs value has been assigned. Conversely, QC flags
  /// previously set to `pass` are reset to `missing` at locations where the variable is set to a
  /// missing value.
  ///
  /// \param varName Variable name.
  /// \param obsVactor Variable values.
  template <typename T>
  void putObservation(const std::string &varName, const std::vector<T> &obsVector,
                      const std::string &outputTag = "DerivedObsValue") {
    std::vector<T> outputObsVector(obsVector);
    // Fill in missing values with values of non-derived group
    if (options_.FillMissingDerivedFromOriginal &&
        outputTag.substr(0, 7) == "Derived") {
      std::string originalTag = outputTag;
      originalTag.erase(0, 7);  // Remove Derived from group name
      if (obsdb_.has(originalTag, varName)) {
        std::vector <T> originalValues(obsdb_.nlocs());
        obsdb_.get_db(originalTag, varName, originalValues);
        const T missing = util::missingValue<T>();
        for (size_t jloc = 0; jloc < obsdb_.nlocs(); ++jloc) {
          if (outputObsVector[jloc] == missing &&
              originalValues[jloc] != missing) {
            outputObsVector[jloc] = originalValues[jloc];
          }
        }
      }
    }
    obsdb_.put_db(outputTag, varName, outputObsVector);

    // Update QC flags to account for values that were previously missing but
    // are now present (or vice versa).
    if (flags_.has(varName)) {
      std::vector<int> &varFlags = flags_[varName];
      ASSERT(varFlags.size() == outputObsVector.size());

      const T missing = util::missingValue<T>();
      for (size_t iloc = 0; iloc < outputObsVector.size(); ++iloc) {
        if (varFlags[iloc] == QCflags::missing && outputObsVector[iloc] != missing)
          varFlags[iloc] = QCflags::pass;
        else if (varFlags[iloc] == QCflags::pass && outputObsVector[iloc] == missing)
          varFlags[iloc] = QCflags::missing;
      }
    }
  }

  /// \brief Save a transformed 2d variable to the `DerivedObsValue` group of the obs space.
  ///
  /// If the saved variable is a simulated variable, QC flags are only set to `pass` at locations
  /// where at least one valid obs value has been assigned.
  ///
  /// \param varName Variable name.
  /// \param channel channel number string.
  /// \param obsVactor Variable values.
  /// \param dimList List of dimension names.
  template <typename T>
  void putObservation(const std::string &varName,
                      const std::string &channel,
                      const std::vector<T> &obsVector,
                      const std::vector<std::string> & dimList,
                      const std::string &outputTag = "DerivedObsValue") {
    std::vector<T> outputObsVector(obsVector);
    // Fill in missing values with values of non-derived group
    if (options_.FillMissingDerivedFromOriginal &&
        outputTag.substr(0, 7) == "Derived") {
      std::string originalTag = outputTag;
      originalTag.erase(0, 7);  // Remove Derived from group name
      if (obsdb_.has(originalTag, varName + "_" + channel)) {
        std::vector <T> originalValues(obsdb_.nlocs());
        // Need to skipDerived so we pickup original ObsValues
        obsdb_.get_db(originalTag, varName + "_" + channel, originalValues, {}, true);
        const T missing = util::missingValue<T>();
        for (size_t jloc = 0; jloc < obsdb_.nlocs(); ++jloc) {
          if (outputObsVector[jloc] == missing &&
              originalValues[jloc] != missing) {
            outputObsVector[jloc] = originalValues[jloc];
          }
        }
      }
    }
    obsdb_.put_db(outputTag, varName + "_" + channel, outputObsVector, dimList);
    if (flags_.has(varName + "_" + channel)) {
      std::vector<int> &varFlags = flags_[varName + "_" + channel];
      ASSERT(varFlags.size() == outputObsVector.size());

      const T missing = util::missingValue<T>();
      for (size_t iloc = 0; iloc < outputObsVector.size(); ++iloc) {
        if (varFlags[iloc] == QCflags::missing && outputObsVector[iloc] != missing)
          varFlags[iloc] = QCflags::pass;
        else if (varFlags[iloc] == QCflags::pass && outputObsVector[iloc] == missing)
          varFlags[iloc] = QCflags::missing;
      }
    }
  }

  std::string getDerivedGroup(const std::string group) const;

  /// subclasses to access Method used for the calculation
  formulas::Method method() const { return method_; }
  bool UseValidDataOnly() const { return UseValidDataOnly_; }
  void SetUseValidDataOnly(bool t) {UseValidDataOnly_ = t; }
  /// subclasses to access the observation name
  std::string obsName() const { return obsName_; }
  /// Configurable parameters
  const VariableTransformParametersBase &options_;

  /// Observation and geoval data
  ObsFilterData data_;
  /// Observation space
  ioda::ObsSpace &obsdb_ = data_.obsspace();
  ioda::ObsDataVector<int> &flags_;
  ioda::ObsDataVector<float> &obserr_;
  /// Missing value (int)
  const int missingValueInt = util::missingValue<int>();
  /// Missing value (float)
  const float missingValueFloat = util::missingValue<float>();
};

/// \brief Transform factory
class TransformFactory {
 public:
  static std::unique_ptr<TransformBase> create(
      const std::string &, const VariableTransformParametersBase &,
      const ObsFilterData &,
      const std::shared_ptr<ioda::ObsDataVector<int>> &,
      const std::shared_ptr<ioda::ObsDataVector<float>> &);
  static std::unique_ptr<VariableTransformParametersBase> createParameters(
      const std::string &);
  virtual ~TransformFactory() = default;

 protected:
  explicit TransformFactory(const std::string &);

 private:
  virtual std::unique_ptr<TransformBase> make(
      const VariableTransformParametersBase &,
      const ObsFilterData &,
      const std::shared_ptr<ioda::ObsDataVector<int>> &,
      const std::shared_ptr<ioda::ObsDataVector<float>> &) = 0;
  virtual std::unique_ptr<VariableTransformParametersBase> makeParameters() const = 0;
  static std::map<std::string, TransformFactory *> &getMakers() {
    static std::map<std::string, TransformFactory *> makers_;
    return makers_;
  }
};

/// \brief Transform maker
template <class T>
class TransformMaker : public TransformFactory {
 private:
  typedef oops::TParameters_IfAvailableElseFallbackType_t<T, GenericVariableTransformParameters>
    Parameters_;

 public:
  explicit TransformMaker(const std::string &name) : TransformFactory(name) {}

  std::unique_ptr<TransformBase> make(
      const VariableTransformParametersBase &options,
      const ObsFilterData& data,
      const std::shared_ptr<ioda::ObsDataVector<int>> &flags,
      const std::shared_ptr<ioda::ObsDataVector<float>> &obserr) override {
    const auto &stronglyTypedParams = dynamic_cast<const Parameters_&>(options);
    return std::unique_ptr<TransformBase>(new T(stronglyTypedParams, data, flags, obserr));
  }

  std::unique_ptr<VariableTransformParametersBase> makeParameters() const override {
    return std::make_unique<Parameters_>();
  }
};

}  // namespace ufo

#endif  // UFO_VARIABLETRANSFORMS_TRANSFORMBASE_H_
