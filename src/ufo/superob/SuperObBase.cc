/*
 * (C) Crown copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/superob/SuperObBase.h"

#include <cstdint>
#include <map>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>

#include "oops/util/Logger.h"
#include "ufo/filters/Variables.h"

namespace ufo {

SuperObBase::SuperObBase(const SuperObParametersBase & /*params*/,
                         const ObsFilterData & data,
                         const std::vector<bool> & apply,
                         const Variables & filtervars,
                         const ioda::ObsDataVector<int> & flags,
                         std::vector<std::vector<bool>> & flagged)
  : data_(data),
    obsdb_(data.obsspace()),
    apply_(apply),
    filtervars_(filtervars),
    flags_(flags),
    flagged_(flagged)
{}

ValidatedSuperObParameters SuperObBase::validateAndParseParameters(
    const SuperObParameters & params) const {
  // Also don't allow variables to increment to be in group ObsValue or
  // DerivedObsValue, since qc flag information cannot be updated.
  if (params.variablesToIncrement.value()) {
    for (const Variable & var : *params.variablesToIncrement.value()) {
      if (var.group() == "ObsValue" || var.group() == "DerivedObsValue") {
        throw eckit::BadParameter(
            "Variables to increment cannot be in group ObsValue or "
            "DerivedObsValue, since qc flag information cannot be updated. "
            "Variable " + var.variable() + " is in group " + var.group(),
            Here());
      }
    }
  }

  // Parse and check parameters
  ValidatedSuperObParameters validated;

  validated.setValuesOutsideWhereClauseToMissing =
      params.setValuesOutsideWhereClauseToMissing.value()
          ? *params.setValuesOutsideWhereClauseToMissing.value()
          : std::vector<bool>(filtervars_.nvars(), true);

  validated.incrementIfNonMissing =
      params.incrementIfNonMissing.value()
          ? *params.incrementIfNonMissing.value()
          : std::vector<bool>(filtervars_.nvars(), false);

  validated.variablesToIncrement =
      params.variablesToIncrement.value()
          ? Variables{*params.variablesToIncrement.value()}
          : Variables{};

  validated.incrementValues = params.incrementValues.value()
                                  ? *params.incrementValues.value()
                                  : std::vector<int>{};

  validated.incrementWholeRecord =
      params.incrementWholeRecord.value()
          ? *params.incrementWholeRecord.value()
          : std::vector<bool>(filtervars_.nvars(), false);

  validated.incrementWholeRecordRespectsWhere =
      params.incrementWholeRecordRespectsWhere.value()
          ? *params.incrementWholeRecordRespectsWhere.value()
          : std::vector<bool>(filtervars_.nvars(), true);

  if (validated.setValuesOutsideWhereClauseToMissing.size() !=
      filtervars_.nvars()) {
    throw eckit::BadParameter(
        "Size of \"set values outside where clause to missing\" does not "
        "match number of filter variables",
        Here());
  }
  if (validated.incrementIfNonMissing.size() != filtervars_.nvars()) {
    throw eckit::BadParameter(
        "Size of \"increment if non-missing\" does not match number of "
        "filter variables",
        Here());
  }
  if (params.incrementIfNonMissing.value()) {
    if (validated.incrementIfNonMissing.size() > 0 &&
        validated.variablesToIncrement.nvars() == 0) {
      throw eckit::BadParameter(
          "Variables to increment must be specified if \"increment if "
          "non-missing\" is used",
          Here());
    }
    if (validated.variablesToIncrement.nvars() > 0 &&
        validated.incrementValues.size() == 0) {
      throw eckit::BadParameter(
          "Increment values must be specified if \"variables to increment\" "
          "are used",
          Here());
    }
    if (validated.variablesToIncrement.nvars() != filtervars_.nvars()) {
      throw eckit::BadParameter(
          "Number of \"variables to increment\" does not match number of "
          "filter variables",
          Here());
    }
    if (validated.incrementValues.size() != filtervars_.nvars()) {
      throw eckit::BadParameter(
          "Number of \"increment values\" does not match number of filter "
          "variables",
          Here());
    }
    if (validated.incrementWholeRecord.size() != filtervars_.nvars()) {
      throw eckit::BadParameter(
          "Size of \"increment whole record\" does not match number of "
          "filter variables",
          Here());
    }
    if (validated.incrementWholeRecordRespectsWhere.size() !=
        filtervars_.nvars()) {
      throw eckit::BadParameter(
          "Size of \"increment whole record respects where\" does not match "
          "number of filter variables",
          Here());
    }
  }
  if (validated.variablesToIncrement.nvars() !=
      validated.incrementValues.size()) {
    throw eckit::BadParameter(
        "Number of \"variables to increment\" does not match number of "
        "\"increment values\"",
        Here());
  }

  // Validate constraints that depend on both parameters and data.
  // For each filter variable, check that referenced variables exist
  // in ObsSpace and that logical constraints are satisfied.
  for (size_t jvar = 0; jvar < filtervars_.nvars(); ++jvar) {
    if (validated.incrementIfNonMissing[jvar]) {
      // Check that the variable to increment exists in ObsSpace
      if (!obsdb_.has(validated.variablesToIncrement[jvar].group(),
                      validated.variablesToIncrement[jvar].variable())) {
        throw eckit::BadParameter(
            "Variable to increment " +
                validated.variablesToIncrement[jvar].variable() + " in group " +
                validated.variablesToIncrement[jvar].group() +
                " does not exist in ObsSpace",
            Here());
      }
      // Check logical constraint: incrementWholeRecordRespectsWhere can only be
      // false if incrementWholeRecord is true
      if (!validated.incrementWholeRecord[jvar] &&
          !validated.incrementWholeRecordRespectsWhere[jvar]) {
        throw eckit::BadParameter(
            "\"increment whole record respects where\" can only be false if "
            "\"increment whole record\" is true",
            Here());
      }
    }
  }

  return validated;
}

void SuperObBase::runAlgorithm(const SuperObParameters & options) const {
  oops::Log::trace() << "SuperObBase::runAlgorithm starting" << std::endl;
  const float missing = util::missingValue<float>();
  const std::size_t nlocs = obsdb_.nlocs();
  const std::vector<std::size_t> & recnums = obsdb_.recidx_all_recnums();

  const ValidatedSuperObParameters validated = validateAndParseParameters(options);

  // Produce lists of locations for each record.
  std::unordered_map<int, std::vector<std::size_t>> locsToUse;

  // Loop over records and fill valid locations in each.
  for (const auto & recnum : recnums) {
    locsToUse[recnum] = {};
    const std::vector<std::size_t> allLocsInRec = obsdb_.recidx_vector(recnum);
    for (size_t jloc : allLocsInRec) {
      if (apply_[jloc]) {
        locsToUse[recnum].push_back(jloc);
        for (size_t jvar = 0; jvar < filtervars_.nvars(); ++jvar) {
          // Set relevant entries in `flagged_` to true. Locations in each
          // record set to `false` indicate where the superob is located.
          flagged_[jvar][jloc] = true;
        }
      }
    }
  }

  // Names of H(x) variables. These are specified regardless of whether
  // H(x) is used in the superob algorithm.
  Variables varhofx(filtervars_, "HofX");

  // Loop over each filter variable and compute superobs for each one.
  // Also compute the associated superob errors, and potentially save any
  // auxiliary variables in the algorithm.
  for (size_t jvar = 0; jvar < filtervars_.nvars(); ++jvar) {
    const std::string variableName = filtervars_[jvar].variable();
    const std::string inputGroupName = filtervars_[jvar].group().empty()
                                           ? "ObsValue"
                                           : filtervars_[jvar].group();
    std::vector<float> obs(nlocs);
    obsdb_.get_db(inputGroupName, variableName, obs);

    std::vector<float> superobs(nlocs, missing);
    if (!validated.setValuesOutsideWhereClauseToMissing[jvar]) {
      superobs = obs;  // preserve obs outside the where clause
      for (const auto& recnum : recnums) {
        for (size_t jloc : locsToUse[recnum]) {
          // still need locs inside where-clause to start as missing
          superobs[jloc] = missing;
        }
      }
    }

    // Vector of H(x) values.
    // Set to missing by default. If `requireHofX()` returns `true`
    // then this vector is filled with H(x) values.
    std::vector<float> hofx(nlocs, missing);
    if (requireHofX()) {
      data_.get(varhofx.variable(jvar), hofx);
    }
    std::vector<int> variableToIncrementData(nlocs, 0);
    if (validated.incrementIfNonMissing[jvar]) {
      obsdb_.get_db(validated.variablesToIncrement[jvar].group(),
                    validated.variablesToIncrement[jvar].variable(),
                    variableToIncrementData);
    }
    for (const auto& recnum : recnums) {
      const bool superObComputed =
          computeSuperOb(locsToUse[recnum], obs, hofx, flags_[variableName],
                         superobs, flagged_[jvar]);
      if (superObComputed && validated.incrementIfNonMissing[jvar]) {
        // Where increment whole record is true, increment all locations
        // in the record selected by the where clause. Otherwise only
        // increment the location(s) where the superob is stored (i.e. where
        // flagged_[jvar] is false).
        const std::vector<std::size_t>& incrementLocs =
            validated.incrementWholeRecordRespectsWhere[jvar]
                ? locsToUse[recnum]  // only locations selected by where clause
                : obsdb_.recidx_vector(
                      recnum);  // all locations in the record (only possible if
                                // incrementWholeRecord[jvar] is true)
        for (size_t jloc : incrementLocs) {
          if (validated.incrementWholeRecord[jvar] || !flagged_[jvar][jloc]) {
            variableToIncrementData[jloc] += validated.incrementValues[jvar];
          }
        }
      }
    }
    // Save the superob values to the ObsSpace.
    obsdb_.put_db("DerivedObsValue", variableName, superobs, filtervars_[jvar].dimList());
    if (validated.incrementIfNonMissing[jvar]) {
      // Save the incremented variable back to the ObsSpace.
      obsdb_.put_db(validated.variablesToIncrement[jvar].group(),
                    validated.variablesToIncrement[jvar].variable(),
                    variableToIncrementData,
                    filtervars_[jvar].dimList());
    }
    // Save any auxiliary variables to the ObsSpace.
    saveAuxiliaryVariables(variableName);
  }
  oops::Log::trace() << "SuperObBase::runAlgorithm done" << std::endl;
}


SuperObFactory::SuperObFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end())
    throw eckit::BadParameter(name + " already registered in ufo::SuperObFactory.", Here());
  getMakers()[name] = this;
}

std::unique_ptr<SuperObBase>
SuperObFactory::create(const SuperObParametersBase & params,
                       const ObsFilterData & data,
                       const std::vector<bool> & apply,
                       const Variables & filtervars,
                       const ioda::ObsDataVector<int> & flags,
                       std::vector<std::vector<bool>> & flagged) {
  oops::Log::trace() << "SuperObBase::create starting" << std::endl;
  const std::string & name = params.superObName;
  typename std::map<std::string, SuperObFactory*>::iterator jloc = getMakers().find(name);
  if (jloc == getMakers().end()) {
    std::string makerNameList;
    for (const auto & makerDetails : getMakers()) makerNameList += "\n  " + makerDetails.first;
    throw eckit::BadParameter(name + " does not exist in ufo::SuperObFactory. "
                              "Possible values:" + makerNameList, Here());
  }
  std::unique_ptr<SuperObBase> ptr =
    jloc->second->make(params, data, apply, filtervars, flags, flagged);
  oops::Log::trace() << "SuperObBase::create done" << std::endl;
  return ptr;
}

std::unique_ptr<SuperObParametersBase>
SuperObFactory::createParameters(const std::string & name) {
  oops::Log::trace() << "SuperObBase::createParameters starting" << std::endl;
  typename std::map<std::string, SuperObFactory*>::iterator jloc = getMakers().find(name);
  if (jloc == getMakers().end()) {
    std::string makerNameList;
    for (const auto & makerDetails : getMakers()) makerNameList += "\n  " + makerDetails.first;
    throw eckit::BadParameter(name + " does not exist in ufo::SuperObFactory. "
                              "Possible values:" + makerNameList, Here());
  }
  std::unique_ptr<SuperObParametersBase> ptr = jloc->second->makeParameters();
  oops::Log::trace() << "SuperObBase::createParameters done" << std::endl;
  return ptr;
}

std::vector<std::size_t> SuperObBase::getUniqueLocations(
    const std::vector<std::size_t>& locs, const Variable& groupingVariable,
    const std::vector<float>& obs, const ioda::ObsDataRow<int>& flags) const {
  switch (obsdb_.dtype(groupingVariable.group(), groupingVariable.variable())) {
    case ioda::ObsDtype::Integer:
      return deduplicateLocations<int>(locs, groupingVariable, obs, flags);
    case ioda::ObsDtype::Integer_64:
      return deduplicateLocations<int64_t>(locs, groupingVariable, obs, flags);
    case ioda::ObsDtype::Float:
      return deduplicateLocations<float>(locs, groupingVariable, obs, flags);
    case ioda::ObsDtype::String:
      return deduplicateLocations<std::string>(locs, groupingVariable, obs,
                                               flags);
    case ioda::ObsDtype::DateTime:
      return deduplicateLocations<util::DateTime>(locs, groupingVariable, obs,
                                                  flags);
    case ioda::ObsDtype::Bool:
      return deduplicateLocations<bool>(locs, groupingVariable, obs, flags);
    default:
      throw eckit::BadParameter(
          "Unsupported grouping variable data type for deduplication", Here());
  }
}

template <typename T>
std::vector<std::size_t> SuperObBase::deduplicateLocations(
    const std::vector<std::size_t>& locs, const Variable& groupingVariable,
    const std::vector<float>& obs, const ioda::ObsDataRow<int>& flags) const {
  // Local hash function for util::DateTime
  // todo: this is a hack - there should be a global hash for DateTime
  struct DateTimeHash {
    std::size_t operator()(const util::DateTime& dt) const {
      return std::hash<std::string>{}(dt.toString());
    }
  };
  using HashType =
      typename std::conditional<std::is_same<T, util::DateTime>::value,
                                DateTimeHash, std::hash<T>>::type;

  std::unordered_map<T, std::pair<float, int>, HashType>
      groupingValueToObsAndFlags;
  std::vector<std::size_t> uniqueLocs;

  std::vector<T> groupingValues;
  obsdb_.get_db(groupingVariable.group(), groupingVariable.variable(),
                groupingValues);

  for (std::size_t jloc : locs) {
    const T& groupingValue = groupingValues[jloc];
    if (groupingValue == util::missingValue<T>()) {
      continue;  // Ignore missing values
    }
    const auto it = groupingValueToObsAndFlags.find(groupingValue);
    if (it == groupingValueToObsAndFlags.end()) {
      groupingValueToObsAndFlags.emplace(
          groupingValue, std::make_pair(obs[jloc], flags[jloc]));
      uniqueLocs.push_back(jloc);
    } else if (obs[jloc] != it->second.first ||
               flags[jloc] != it->second.second) {
      throw eckit::UserError(
          "Grouping variable has non-identical observation values or QC flags "
          "within a group.",
          Here());
    }
  }

  return uniqueLocs;
}

void SuperObBase::assignSuperObToLocations(const std::vector<std::size_t>& locs,
                                           const float superobValue,
                                           const std::size_t superobloc,
                                           const bool assignToAll,
                                           std::vector<float>& superobs,
                                           std::vector<bool>& flagged) const {
  const float missing = util::missingValue<float>();

  for (std::size_t jloc : locs) {
    if (assignToAll || jloc == superobloc) {
      superobs[jloc] = superobValue;
      flagged[jloc] = false;
    } else {
      superobs[jloc] = missing;
      flagged[jloc] = true;
    }
  }
}

// Explicit instantiations for the supported grouping-variable types
template std::vector<std::size_t> SuperObBase::deduplicateLocations<int>(
    const std::vector<std::size_t>&, const Variable&, const std::vector<float>&,
    const ioda::ObsDataRow<int>&) const;

template std::vector<std::size_t> SuperObBase::deduplicateLocations<int64_t>(
    const std::vector<std::size_t>&, const Variable&, const std::vector<float>&,
    const ioda::ObsDataRow<int>&) const;

template std::vector<std::size_t> SuperObBase::deduplicateLocations<float>(
    const std::vector<std::size_t>&, const Variable&, const std::vector<float>&,
    const ioda::ObsDataRow<int>&) const;

template std::vector<std::size_t> SuperObBase::deduplicateLocations<
    std::string>(const std::vector<std::size_t>&, const Variable&,
                 const std::vector<float>&, const ioda::ObsDataRow<int>&) const;

template std::vector<std::size_t>
SuperObBase::deduplicateLocations<util::DateTime>(
    const std::vector<std::size_t>&, const Variable&, const std::vector<float>&,
    const ioda::ObsDataRow<int>&) const;

template std::vector<std::size_t> SuperObBase::deduplicateLocations<bool>(
    const std::vector<std::size_t>&, const Variable&, const std::vector<float>&,
    const ioda::ObsDataRow<int>&) const;
}  // namespace ufo
