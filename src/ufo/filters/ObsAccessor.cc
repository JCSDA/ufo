/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/ObsAccessor.h"

#include <memory>
#include <string>
#include <vector>

#include "ioda/distribution/InefficientDistribution.h"
#include "ioda/ObsSpace.h"
#include "oops/mpi/mpi.h"

#include "ufo/filters/FilterUtils.h"
#include "ufo/filters/QCflags.h"
#include "ufo/filters/Variables.h"
#include "ufo/utils/RecursiveSplitter.h"

namespace ufo {

namespace {

template <typename VariableType>
std::vector<VariableType> getOnlySelectedVariable(const std::string &group,
                                                  const std::string &variable,
                                                  const std::vector<bool> &apply,
                                                  const std::vector<size_t> &glocs,
                                                  const ioda::ObsSpace &obsdb,
                                                  const ioda::Distribution &obsDistribution) {
  const size_t nlocs = obsdb.nlocs();
  const size_t gnlocs = obsdb.globalNumLocs();

  // Local result vector.
  std::vector<VariableType> resultLocal(nlocs);
  obsdb.get_db(group, variable, resultLocal);

  // Vector signifying if each location corresponds to a 'patch' observation.
  // This is relevant when locations are held on more than one rank, because
  // it enables a single rank to be assigned to each location.
  std::vector<bool> patchObsVec(nlocs);
  obsDistribution.patchObs(patchObsVec);

  // Make local result vector smaller according to both `apply` and `patchObsVec`.
  std::vector<VariableType> resultReduced;
  for (size_t i = 0; i < nlocs; ++i) {
    if (apply[i] && patchObsVec[i]) {
      resultReduced.push_back(resultLocal[i]);
    }
  }

  // Gather reduced result vector.
  // Note: `obsDistribution.allGatherv()` throws an assertion if `resultReduced.size() != nlocs`,
  // which will occur if at least one entry in `apply` is `false`.
  // Therefore use the underlying `oops::mpi::allGatherv` instead.
  oops::mpi::allGatherv(obsdb.comm(), resultReduced);
  ASSERT(resultReduced.size() == glocs.size());

  // Set up a global result vector that is initially filled with missing values.
  const VariableType missing = util::missingValue<VariableType>();
  std::vector<VariableType> result(gnlocs, missing);
  // Put the gathered result into the global result vector according to the `glocs` vector.
  for (size_t i = 0; i < resultReduced.size(); ++i) {
    result[glocs[i]] = resultReduced[i];
  }

  return result;
}

template <typename VariableType>
std::vector<VariableType> getVariableFromObsSpaceImpl(
    const std::string &group, const std::string &variable,
    const ioda::ObsSpace &obsdb, const ioda::Distribution &obsDistribution) {
  std::vector<VariableType> result(obsdb.nlocs());
  obsdb.get_db(group, variable, result);
  obsDistribution.allGatherv(result);
  return result;
}

template <typename VariableType>
std::vector<VariableType> getVariableFromObsSpaceImpl(
    const std::string &group,
    const std::string &variable,
    const std::vector<bool> &apply,
    const std::vector<size_t> &glocs,
    const bool groupByCategoryVariable,
    const ioda::ObsSpace &obsdb,
    const ioda::Distribution &obsDistribution,
    const bool accountForWhere) {

  // If each record is held by a single rank, no MPI communication is performed
  // so there is no need to gather observations across ranks.
  if (groupByCategoryVariable || !accountForWhere) {
    return getVariableFromObsSpaceImpl<VariableType>(group, variable, obsdb, obsDistribution);
  }

  return getOnlySelectedVariable<VariableType>(group,
                                               variable,
                                               apply,
                                               glocs,
                                               obsdb,
                                               obsDistribution);
}

/// Return the vector of elements of \p categories with indices \p validObsIds.
template <typename T>
std::vector<T> getValidObservationCategories(const std::vector<T> &categories,
                                             const std::vector<size_t> &validObsIds) {
  std::vector<T> validObsCategories(validObsIds.size());
  for (size_t validObsIndex = 0; validObsIndex < validObsIds.size(); ++validObsIndex) {
    validObsCategories[validObsIndex] = categories[validObsIds[validObsIndex]];
  }
  return validObsCategories;
}

template <typename VariableType>
void groupObservationsByVariableImpl(
    const Variable &variable,
    const std::vector<size_t> &validObsIds,
    const ioda::ObsSpace &obsdb,
    const ioda::Distribution &obsDistribution,
    RecursiveSplitter &splitter) {
  std::vector<VariableType> obsCategories(obsdb.nlocs());
  obsdb.get_db(variable.group(), variable.variable(), obsCategories);
  obsDistribution.allGatherv(obsCategories);
  const std::vector<VariableType> validObsCategories = getValidObservationCategories(
        obsCategories, validObsIds);

  splitter.groupBy(validObsCategories);
}

template <typename VariableType>
void groupObservationsByVariableImpl(
    const Variable &variable,
    const std::vector<bool> &apply,
    const std::vector<size_t> &glocs,
    const std::vector<size_t> &validObsIds,
    const ioda::ObsSpace &obsdb,
    const ioda::Distribution &obsDistribution,
    RecursiveSplitter &splitter,
    const bool accountForWhere) {

  std::vector<VariableType> obsCategories;
  if (!accountForWhere) {
    obsdb.get_db(variable.group(), variable.variable(), obsCategories);
    obsDistribution.allGatherv(obsCategories);
  } else {
    obsCategories = getOnlySelectedVariable<VariableType>(variable.group(),
                                                          variable.variable(),
                                                          apply,
                                                          glocs,
                                                          obsdb,
                                                          obsDistribution);
  }

  const std::vector<VariableType> validObsCategories = getValidObservationCategories(
        obsCategories, validObsIds);

  splitter.groupBy(validObsCategories);
}

}  // namespace

ObsAccessor::ObsAccessor(const ioda::ObsSpace &obsdb,
                         GroupBy groupBy,
                         boost::optional<Variable> categoryVariable,
                         const bool accountForWhere)
  : obsdb_(&obsdb), groupBy_(groupBy), categoryVariable_(categoryVariable),
    accountForWhere_(accountForWhere)
{
  oops::Log::trace() << "ObservationAccessor constructor" << std::endl;
  // If the observations are to be grouped by a category variable, and that variable was
  // also used to divide the ObsSpace into records, change the value of `groupBy_`.
  // Do the same if observations are to be grouped by the ObsSpace record index.
  // This is not done if the records are treated as single observations (for which
  // `groupBy_` is equal to `GroupBy::SINGLE_OBS`).
  if ((groupBy_ == GroupBy::VARIABLE && wereRecordsGroupedByCategoryVariable()) ||
      groupBy_ == GroupBy::RECORD_ID) {
    groupBy_ = GroupBy::CATEGORY_VARIABLE;
  }

  if (groupBy_ == GroupBy::CATEGORY_VARIABLE) {
    // Each record is held by a single process, so there's no need to exchange data between
    // processes and we can use an InefficientDistribution rather than the distribution taken from
    // obsdb_. Which in this case is *efficient*!
    obsDistribution_ = std::make_shared<ioda::InefficientDistribution>(obsdb_->comm(),
                                                        ioda::EmptyDistributionParameters());
    oops::Log::trace() << "ObservationAccessor: no MPI communication necessary" << std::endl;
  } else {
    obsDistribution_ = obsdb.distribution();
  }
}

ObsAccessor ObsAccessor::toAllObservations
(const ioda::ObsSpace &obsdb, const bool accountForWhere) {
  return ObsAccessor(obsdb, GroupBy::NOTHING, boost::none, accountForWhere);
}

ObsAccessor ObsAccessor::toObservationsSplitIntoIndependentGroupsByRecordId(
  const ioda::ObsSpace &obsdb, const bool accountForWhere) {
  return ObsAccessor(obsdb, GroupBy::CATEGORY_VARIABLE, boost::none, accountForWhere);
}

ObsAccessor ObsAccessor::toObservationsSplitIntoIndependentGroupsByVariable(
  const ioda::ObsSpace &obsdb, const Variable &variable, const bool accountForWhere) {
  return ObsAccessor(obsdb, GroupBy::VARIABLE, variable, accountForWhere);
}

ObsAccessor ObsAccessor::toSingleObservationsSplitIntoIndependentGroupsByVariable(
  const ioda::ObsSpace &obsdb, const Variable &variable, const bool accountForWhere) {
  return ObsAccessor(obsdb, GroupBy::SINGLE_OBS, variable, accountForWhere);
}

std::vector<bool> ObsAccessor::getGlobalApply(
    const std::vector<bool> &apply) const {
  std::vector<int> globalApply(apply.begin(), apply.end());
  obsDistribution_->allGatherv(globalApply);
  return std::vector<bool>(globalApply.begin(), globalApply.end());
}

std::vector<size_t> ObsAccessor::getValidObservationIds(
    const std::vector<bool> &apply, const ioda::ObsDataVector<int> &flags,
    const ufo::Variables &filtervars, bool candidateForRetentionIfAnyFilterVariablesPassedQC)
    const {
  std::vector<bool> isValid = apply;
  const UnselectLocationIf mode = candidateForRetentionIfAnyFilterVariablesPassedQC ?
        UnselectLocationIf::ALL_FILTER_VARIABLES_REJECTED :
        UnselectLocationIf::ANY_FILTER_VARIABLE_REJECTED;
  unselectRejectedLocations(isValid, filtervars, flags, mode);

  std::vector<int> globalIsValid;

  if (groupBy_ == GroupBy::CATEGORY_VARIABLE ||
      !accountForWhere_) {
    // If observations are grouped according to the category variable
    // can simply gather the globalIsValid vector using the underlying ObsSpace
    // distribution.
    // TODO(wsmigaj): use std::vector<unsigned char> to save space
    globalIsValid.insert(globalIsValid.begin(), isValid.cbegin(), isValid.cend());
    obsDistribution_->allGatherv(globalIsValid);
  } else {
    const size_t nlocs = obsdb_->nlocs();
    const size_t gnlocs = obsdb_->globalNumLocs();
    // Vector signifying if each location corresponds to a 'patch' observation.
    // This is relevant when locations are held on more than one rank, because
    // it enables a single rank to be assigned to each location.
    std::vector<bool> patchObsVec(nlocs);
    obsDistribution_->patchObs(patchObsVec);
    std::vector<size_t> glocs;
    for (size_t i = 0; i < nlocs; ++i) {
      if (isValid[i] && patchObsVec[i]) {
        const size_t gloc = obsDistribution_->globalUniqueConsecutiveLocationIndex(i);
        glocs.push_back(gloc);
      }
    }
    // Gather `glocs` vector. If any entries in `isValid` are `false`, the length
    // of `glocs` will be smaller than the global number of locations.
    oops::mpi::allGatherv(obsdb_->comm(), glocs);
    globalIsValid.assign(gnlocs, 0);
    for (size_t i = 0; i < glocs.size(); ++i) {
      globalIsValid[glocs[i]] = 1;
    }
  }

  std::vector<size_t> validObsIds;
  for (size_t obsId = 0; obsId < globalIsValid.size(); ++obsId)
    if (globalIsValid[obsId])
      validObsIds.push_back(obsId);

  return validObsIds;
}

std::vector<size_t> ObsAccessor::getValidObservationIds(
    const std::vector<bool> &apply) const {
  // TODO(wsmigaj): use std::vector<unsigned char> to save space
  std::vector<int> globalIsValid(apply.begin(), apply.end());
  obsDistribution_->allGatherv(globalIsValid);

  std::vector<size_t> validObsIds;
  for (size_t obsId = 0; obsId < globalIsValid.size(); ++obsId)
    if (globalIsValid[obsId])
      validObsIds.push_back(obsId);

  return validObsIds;
}


/// Get valid (non-missing, where-included) obs indices for a given profile.
const std::vector<size_t> ObsAccessor::getValidObsIdsInProfile(const size_t & iProfile,
                                      const std::vector<bool> & apply,
                                      const ioda::ObsDataVector<int> &flags,
                                      const Variables & filtervars,
                                      bool candidateForRetentionIfAnyFilterVariablesPassedQC)
                                      const {
  // Get vector of obs within the profile:
  const std::vector<size_t> & obs_inds = obsdb_->recidx_vector(iProfile);
  std::vector<size_t> validObsIds;
  std::vector<bool> isValid = apply;
  const UnselectLocationIf mode = candidateForRetentionIfAnyFilterVariablesPassedQC ?
        UnselectLocationIf::ALL_FILTER_VARIABLES_REJECTED :
        UnselectLocationIf::ANY_FILTER_VARIABLE_REJECTED;
  unselectRejectedLocations(isValid, filtervars, flags, mode, obs_inds);
  for (size_t ind = 0; ind < obs_inds.size(); ++ind) {
    if (isValid[obs_inds[ind]]) {
      validObsIds.push_back(obs_inds[ind]);
    }
  }
  return validObsIds;
}


std::vector<int> ObsAccessor::getIntVariableFromObsSpace(
    const std::string &group, const std::string &variable) const {
  return getVariableFromObsSpaceImpl<int>(group, variable, *obsdb_, *obsDistribution_);
}

std::vector<float> ObsAccessor::getFloatVariableFromObsSpace(
    const std::string &group, const std::string &variable) const {
  return getVariableFromObsSpaceImpl<float>(group, variable, *obsdb_, *obsDistribution_);
}

std::vector<double> ObsAccessor::getDoubleVariableFromObsSpace(
    const std::string &group, const std::string &variable) const {
  return getVariableFromObsSpaceImpl<double>(group, variable, *obsdb_, *obsDistribution_);
}

std::vector<std::string> ObsAccessor::getStringVariableFromObsSpace(
    const std::string &group, const std::string &variable) const {
  return getVariableFromObsSpaceImpl<std::string>(group, variable, *obsdb_, *obsDistribution_);
}

std::vector<util::DateTime> ObsAccessor::getDateTimeVariableFromObsSpace(
      const std::string &group, const std::string &variable) const {
  return getVariableFromObsSpaceImpl<util::DateTime>(group, variable, *obsdb_, *obsDistribution_);
}

std::vector<int> ObsAccessor::getIntVariableFromObsSpace
(const std::string &group,
 const std::string &variable,
 const std::vector<bool> &apply) const {
  const bool groupByCategoryVariable = groupBy_ == GroupBy::CATEGORY_VARIABLE;
  if (!groupByCategoryVariable) {
    fillGlocs(apply);
  }
  return getVariableFromObsSpaceImpl<int>(group, variable, apply, glocs_,
                                          groupByCategoryVariable,
                                          *obsdb_, *obsDistribution_,
                                          accountForWhere_);
}

std::vector<float> ObsAccessor::getFloatVariableFromObsSpace
(const std::string &group,
 const std::string &variable,
 const std::vector<bool> &apply) const {
  const bool groupByCategoryVariable = groupBy_ == GroupBy::CATEGORY_VARIABLE;
  if (!groupByCategoryVariable) {
    fillGlocs(apply);
  }
  return getVariableFromObsSpaceImpl<float>(group, variable, apply, glocs_,
                                            groupByCategoryVariable,
                                            *obsdb_, *obsDistribution_,
                                            accountForWhere_);
}

std::vector<util::DateTime> ObsAccessor::getDateTimeVariableFromObsSpace
(const std::string &group,
 const std::string &variable,
 const std::vector<bool> &apply) const {
  const bool groupByCategoryVariable = groupBy_ == GroupBy::CATEGORY_VARIABLE;
  if (!groupByCategoryVariable) {
    fillGlocs(apply);
  }
  return getVariableFromObsSpaceImpl<util::DateTime>(group, variable, apply, glocs_,
                                                     groupByCategoryVariable,
                                                     *obsdb_, *obsDistribution_,
                                                     accountForWhere_);
}

std::vector<size_t> ObsAccessor::getRecordIds() const {
  std::vector<size_t> recordIds = obsdb_->recnum();
  obsDistribution_->allGatherv(recordIds);
  return recordIds;
}

std::vector<bool> ObsAccessor::getBoolVariableFromObsSpace(
  const std::string &group, const std::string &variable) const {
  std::vector<bool> requestedVariable;
  obsdb_->get_db(group, variable, requestedVariable);
  std::vector<int> globalRequestedVariable(requestedVariable.begin(), requestedVariable.end());
  obsDistribution_->allGatherv(globalRequestedVariable);
  std::vector<bool> globalRequestedVariableBool(globalRequestedVariable.begin(),
                                                globalRequestedVariable.end());
  return globalRequestedVariableBool;
}

size_t ObsAccessor::totalNumObservations() const {
  return obsdb_->globalNumLocs();
}

RecursiveSplitter ObsAccessor::splitObservationsIntoIndependentGroups(
    const std::vector<size_t> &validObsIds, bool opsCompatibilityMode) const {
  RecursiveSplitter splitter(validObsIds.size(), opsCompatibilityMode);
  switch (groupBy_) {
  case GroupBy::NOTHING:
    // Nothing to do
    break;
  case GroupBy::CATEGORY_VARIABLE:
    groupObservationsByRecordNumber(validObsIds, splitter);
    break;
  case GroupBy::VARIABLE:
    groupObservationsByCategoryVariable(validObsIds, splitter);
    break;
  case GroupBy::SINGLE_OBS:
    groupObservationsByCategoryVariable(validObsIds, splitter);
    break;
  }
  return splitter;
}

void ObsAccessor::groupObservationsByRecordNumber(const std::vector<size_t> &validObsIds,
                                                           RecursiveSplitter &splitter) const {
  const std::vector<size_t> &obsCategories = obsdb_->recnum();
  std::vector<size_t> validObsCategories = getValidObservationCategories(
        obsCategories, validObsIds);
  splitter.groupBy(validObsCategories);
}

void ObsAccessor::groupObservationsByCategoryVariable(
    const std::vector<size_t> &validObsIds,
    RecursiveSplitter &splitter) const {
  switch (obsdb_->dtype(categoryVariable_->group(), categoryVariable_->variable())) {
  case ioda::ObsDtype::Integer:
    groupObservationsByVariableImpl<int>(*categoryVariable_, validObsIds,
                                         *obsdb_, *obsDistribution_, splitter);
    break;

  case ioda::ObsDtype::String:
    groupObservationsByVariableImpl<std::string>(*categoryVariable_, validObsIds,
                                                 *obsdb_, *obsDistribution_, splitter);
    break;

  default:
    throw eckit::UserError(
          categoryVariable_->group() + "/" + categoryVariable_->variable() +
          " is neither an integer nor a string variable", Here());
  }
}

void ObsAccessor::flagRejectedObservations(
    const std::vector<bool> &isRejected, std::vector<std::vector<bool> > &flagged) const {
  const size_t localNumObs = obsdb_->nlocs();
  for (const std::vector<bool> & variableFlagged : flagged)
    ASSERT(variableFlagged.size() == localNumObs);

  for (size_t localObsId = 0; localObsId < localNumObs; ++localObsId) {
    const size_t globalObsId =
        obsDistribution_->globalUniqueConsecutiveLocationIndex(localObsId);
    if (isRejected[globalObsId]) {
      for (std::vector<bool> & variableFlagged : flagged)
        variableFlagged[localObsId] = true;
    }
  }
}

void ObsAccessor::flagObservationsForAnyFilterVariableFailingQC(
    const std::vector<bool> &apply, const ioda::ObsDataVector<int> &flags,
    const ufo::Variables &filtervars, std::vector<std::vector<bool> > &flagged) const {
  std::vector<size_t> indexOfFilterVariableInFlags;
  for (size_t ivar = 0; ivar < flagged.size(); ++ivar) {
    std::string filterVariableName = filtervars.variable(ivar).variable();
    indexOfFilterVariableInFlags.push_back(flags.varnames().find(filterVariableName));
  }
  for (size_t iloc = 0; iloc < obsdb_->nlocs(); ++iloc) {
    if (apply[iloc]) {
      bool atLeastOneFilterVariableFailsQC = false;
      for (size_t ivar = 0; ivar < flagged.size(); ++ivar) {
        if (QCflags::isRejected(flags[indexOfFilterVariableInFlags[ivar]][iloc])) {
          atLeastOneFilterVariableFailsQC = true;
          break;
        }
      }
      if (atLeastOneFilterVariableFailsQC) {
        for (size_t ivar = 0; ivar < flagged.size(); ++ivar)
          flagged[ivar][iloc] = true;
      }
    }
  }
}

bool ObsAccessor::wereRecordsGroupedByCategoryVariable() const {
  std::vector<std::string> groupingVars = obsdb_->obs_group_vars();
  std::string groupingVar;
  if (groupingVars.size() > 0) {
    groupingVar = groupingVars[0];
  }
  return categoryVariable_ != boost::none &&
         categoryVariable_->variable() == groupingVar &&
         categoryVariable_->group() == "MetaData";
}

}  // namespace ufo
