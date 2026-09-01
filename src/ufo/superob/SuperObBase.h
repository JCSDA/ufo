/*
 * (C) Crown copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_SUPEROB_SUPEROBBASE_H_
#define UFO_SUPEROB_SUPEROBBASE_H_

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/util/AssociativeContainers.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "oops/util/parameters/HasParameters_.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/filters/FilterParametersBase.h"

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/QCflags.h"
#include "ufo/filters/SuperObParameters.h"
#include "ufo/filters/Variables.h"


namespace ufo {

/// \brief SuperOb parameters base class.
class SuperObParametersBase : public oops::Parameters {
  OOPS_ABSTRACT_PARAMETERS(SuperObParametersBase, oops::Parameters)
 public:
  /// Name of the superobbing algorithm.
  /// Valid names are specified using a `SuperObMaker` in subclasses of SuperObBase in the
  /// ufo/superob directory.
  oops::RequiredParameter<std::string> superObName{"name", this};
};

/// \brief Concrete class containing the options specified by the SuperObParametersBase.
/// This can be overwritten by specific superob algorithms if these options are not
/// appropriate.
class GenericSuperObParameters : public SuperObParametersBase {
  OOPS_CONCRETE_PARAMETERS(GenericSuperObParameters, SuperObParametersBase)
 public:
  /// If true, the superob value is assigned to all values in the record and
  /// no rejection is performed
  oops::Parameter<bool> assignToAllValues{"assign to all values in record", false, this};

  /// Optional variable used to group duplicate observations within each record.
  /// When specified, locations within a record that share the same value in this
  /// variable are treated as duplicates, and only one representative location from
  /// each unique value is used in the superob calculation (but not assignment).
  /// For example, if the variable contains [0, 0, 0, 1, 1] for observation values
  /// [0.1, 0.1, 0.1, 0.2, 0.2], the superob is computed using the deduplicated set
  /// [0.1, 0.2], but if "assign to all values in record" is true, all 5 observations
  /// receive the computed superob value. NOTE: if grouped values are not identical
  /// or have differing QC flags within a group, an exception is thrown.
  oops::OptionalParameter<Variable> groupingVariable{
      "grouping variable", this};
};

/// \brief Validated and parsed SuperOb parameters.
/// Encapsulates the parsed and validated parameters that will be used throughout
/// the superobbing routine.
struct ValidatedSuperObParameters {
  std::vector<bool> setValuesOutsideWhereClauseToMissing;
  std::vector<bool> incrementIfNonMissing;
  Variables variablesToIncrement;
  std::vector<int> incrementValues;
  std::vector<bool> incrementWholeRecord;
  std::vector<bool> incrementWholeRecordRespectsWhere;
};

/// \brief SuperOb base class.
///
/// Subclasses of this class must implement the `computeSuperOb` method.
/// The `requireHofX` function returns `true` by default. If H(x) is not required
/// by a subclass then the implementation of `requireHofX` in that subclass
/// should be set to return `false`.
/// If required, the subclass should also override the `saveAuxiliaryVariables` method,
/// which can be used to save any desired quantities to the ObsSpace.
class SuperObBase {
 public:
  explicit SuperObBase(const SuperObParametersBase &,
                       const ObsFilterData &,
                       const std::vector<bool> &,
                       const Variables &,
                       const ioda::ObsDataVector<int> &,
                       std::vector<std::vector<bool>> &);
  virtual ~SuperObBase() {}

  /// Run the chosen superobbing algorithm on each record in the data set.
  void runAlgorithm(const SuperObParameters &) const;

  /// Specifies whether the chosen superobbing algorithm requires H(x).
  virtual bool requireHofX() const {return true;}

 protected:
  /// Compute superob values and errors for each record. Returns `true` if a
  /// superob was successfully computed for the record, `false` otherwise.
  virtual bool computeSuperOb(const std::vector<std::size_t> &,
                              const std::vector<float> &,
                              const std::vector<float> &,
                              const ioda::ObsDataRow<int> &,
                              std::vector<float> &,
                              std::vector<bool> &) const = 0;

  /// Save any auxiliary variables to the ObsSpace.
  /// By default, no extra variables are saved.
  /// The parameter `variableName` is the name of the filter variable.
  virtual void saveAuxiliaryVariables(const std::string & variableName) const {}

  /// Helper methods to deduplicate locations based on grouping variable values.
  /// Returns a vector of locations containing only the first occurrence of each
  /// unique grouping value.
  std::vector<std::size_t> getUniqueLocations(
      const std::vector<std::size_t>& locs,
      const Variable& groupingVariable,
      const std::vector<float>& obs,
      const ioda::ObsDataRow<int>& flags) const;
  /// Deduplicate locs by groupingValues, returning the index of each
  /// unique groupingValue's representative location. For example, given
  /// locs [0, 1, 2, 3] with groupingValues [A, A, B, B], the result
  /// would be [0, 2] (one loc per distinct groupingValue). For duplicate
  /// members of a group, the values in obs and flags must match;
  /// otherwise an exception is thrown.
  template<typename T>
  std::vector<std::size_t> deduplicateLocations(
      const std::vector<std::size_t>& locs,
      const Variable& groupingValues,
      const std::vector<float>& obs,
      const ioda::ObsDataRow<int>& flags) const;

  /// Helper method to assign superob value to locations in a record.
  /// If assignToAll is true, assigns to all locations and unflags them all.
  /// Otherwise, assigns only to superobloc then flags and sets all other
  /// locations to missing.
  void assignSuperObToLocations(const std::vector<std::size_t> &locs,
                                const float superobValue,
                                const std::size_t superobloc,
                                const bool assignToAll,
                                std::vector<float> &superobs,
                                std::vector<bool> &flagged) const;

  const ObsFilterData & data_;
  ioda::ObsSpace & obsdb_;

 private:
  /// Validate and parse SuperOb parameters into a single struct.
  /// Performs all parameter consistency checks, applies defaults, and validates
  /// constraints that depend on both parameters and data (e.g., variable existence).
  /// Throws eckit::BadParameter if validation fails.
  ValidatedSuperObParameters validateAndParseParameters(
      const SuperObParameters & params) const;

  const std::vector<bool> apply_;
  const Variables & filtervars_;
  const ioda::ObsDataVector<int> & flags_;
  std::vector<std::vector<bool>> & flagged_;
};

/// SuperOb factory.
class SuperObFactory {
 public:
  static std::unique_ptr<SuperObBase> create(const SuperObParametersBase &,
                                             const ObsFilterData &,
                                             const std::vector<bool> &,
                                             const Variables &,
                                             const ioda::ObsDataVector<int> &,
                                             std::vector<std::vector<bool>> &);

  static std::unique_ptr<SuperObParametersBase> createParameters(const std::string &name);

  /// \brief Return the names of all predictors that can be created by one of the registered makers.
  static std::vector<std::string> getMakerNames() {
    return oops::keys(getMakers());
  }

  virtual ~SuperObFactory() = default;

 protected:
  explicit SuperObFactory(const std::string &);

 private:
  virtual std::unique_ptr<SuperObBase> make(const SuperObParametersBase &,
                                            const ObsFilterData &,
                                            const std::vector<bool> &,
                                            const Variables &,
                                            const ioda::ObsDataVector<int> &,
                                            std::vector<std::vector<bool>> &) = 0;

  virtual std::unique_ptr<SuperObParametersBase> makeParameters() const = 0;

  static std::map <std::string, SuperObFactory*> & getMakers() {
    static std::map <std::string, SuperObFactory*> makers_;
    return makers_;
  }
};

template<class T>
class SuperObMaker : public SuperObFactory {
  typedef oops::TParameters_IfAvailableElseFallbackType_t<T, GenericSuperObParameters>
    Parameters_;

  std::unique_ptr<SuperObBase> make(const SuperObParametersBase & params,
                                    const ObsFilterData & data,
                                    const std::vector<bool> & apply,
                                    const Variables & filtervars,
                                    const ioda::ObsDataVector<int> & flags,
                                    std::vector<std::vector<bool>> & flagged) override {
    const auto & stronglyTypedParams = dynamic_cast<const Parameters_&>(params);
    return std::unique_ptr<SuperObBase>
      (new T(stronglyTypedParams, data, apply, filtervars, flags, flagged));
  }

  std::unique_ptr<SuperObParametersBase> makeParameters() const override {
    return std::make_unique<Parameters_>();
  }

 public:
  explicit SuperObMaker(const std::string & name)
    : SuperObFactory(name) {}
};

}  // namespace ufo

#endif  // UFO_SUPEROB_SUPEROBBASE_H_
