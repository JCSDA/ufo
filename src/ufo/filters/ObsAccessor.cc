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

#include <boost/make_unique.hpp>

#include "ioda/distribution/InefficientDistribution.h"
#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "ufo/filters/QCflags.h"
#include "ufo/utils/RecursiveSplitter.h"

namespace ufo {

namespace {

template <typename VariableType>
std::vector<VariableType> getVariableFromObsSpaceImpl(
    const std::string &group, const std::string &variable,
    const ioda::ObsSpace &obsdb, const ioda::Distribution &obsDistribution) {
  std::vector<VariableType> result(obsdb.nlocs());
  obsdb.get_db(group, variable, result);
  obsDistribution.allGatherv(result);
  return result;
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

}  // namespace

ObsAccessor::ObsAccessor(const ioda::ObsSpace &obsdb,
                         GroupBy groupBy,
                         boost::optional<Variable> categoryVariable)
  : obsdb_(&obsdb), groupBy_(groupBy), categoryVariable_(categoryVariable)
{
  if (groupBy_ == GroupBy::VARIABLE && wereRecordsGroupedByCategoryVariable())
    groupBy_ = GroupBy::RECORD_ID;

  if (groupBy_ == GroupBy::RECORD_ID) {
    // Each record is held by a single process, so there's no need to exchange data between
    // processes and we can use an InefficientDistribution rather than the distribution taken from
    // obsdb_. Which in this case is *efficient*!
    inefficientDistribution_ = boost::make_unique<ioda::InefficientDistribution>(obsdb_->comm());
    obsDistribution_ = inefficientDistribution_.get();
    oops::Log::trace() << "ObservationAccessor: no MPI communication necessary" << std::endl;
  } else {
    obsDistribution_ = &obsdb.distribution();
  }
}

ObsAccessor::~ObsAccessor() {
  // Defined here rather than in the header file, since here the definition of
  // InefficientDistribution is visible and so inefficientDistribution_ can be correctly destructed.
}

ObsAccessor ObsAccessor::toAllObservations(
    const ioda::ObsSpace &obsdb) {
  return ObsAccessor(obsdb, GroupBy::NOTHING, boost::none);
}

ObsAccessor ObsAccessor::toObservationsSplitIntoIndependentGroupsByRecordId(
    const ioda::ObsSpace &obsdb) {
  return ObsAccessor(obsdb, GroupBy::RECORD_ID, boost::none);
}

ObsAccessor ObsAccessor::toObservationsSplitIntoIndependentGroupsByVariable(
    const ioda::ObsSpace &obsdb, const Variable &variable) {
  return ObsAccessor(obsdb, GroupBy::VARIABLE, variable);
}

std::vector<size_t> ObsAccessor::getValidObservationIds(
    const std::vector<bool> &apply, const ioda::ObsDataVector<int> &flags) const {
  size_t obsIdDisplacement = obsdb_->nlocs();
  obsDistribution_->exclusiveScan(obsIdDisplacement);

  std::vector<size_t> validObsIds;
  for (size_t obsId = 0; obsId < apply.size(); ++obsId)
    if (apply[obsId] && flags[0][obsId] == QCflags::pass)
      validObsIds.push_back(obsIdDisplacement + obsId);

  obsDistribution_->allGatherv(validObsIds);

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

std::vector<size_t> ObsAccessor::getRecordIds() const {
  std::vector<size_t> recordIds = obsdb_->recnum();
  obsDistribution_->allGatherv(recordIds);
  return recordIds;
}

size_t ObsAccessor::totalNumObservations() const {
  size_t totalNumObs = obsdb_->nlocs();
  obsDistribution_->sum(totalNumObs);
  return totalNumObs;
}

RecursiveSplitter ObsAccessor::splitObservationsIntoIndependentGroups(
    const std::vector<size_t> &validObsIds) const {
  RecursiveSplitter splitter(validObsIds.size());
  switch (groupBy_) {
  case GroupBy::NOTHING:
    // Nothing to do
    break;
  case GroupBy::RECORD_ID:
    groupObservationsByRecordNumber(validObsIds, splitter);
    break;
  case GroupBy::VARIABLE:
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
          categoryVariable_->variable() + "@" + categoryVariable_->group() +
          " is neither an integer nor a string variable", Here());
  }
}

void ObsAccessor::flagRejectedObservations(
    const std::vector<bool> &isRejected, std::vector<std::vector<bool> > &flagged) const {
  const size_t localNumObs = obsdb_->nlocs();
  size_t displacement = localNumObs;
  obsDistribution_->exclusiveScan(displacement);

  for (std::vector<bool> & variableFlagged : flagged) {
    ASSERT(variableFlagged.size() == localNumObs);
    for (size_t localObsId = 0; localObsId < localNumObs; ++localObsId) {
      const size_t globalObsId = displacement + localObsId;
      if (isRejected[globalObsId])
        variableFlagged[localObsId] = true;
    }
  }
}

bool ObsAccessor::wereRecordsGroupedByCategoryVariable() const {
  return categoryVariable_ != boost::none &&
         categoryVariable_->variable() == obsdb_->obs_group_var() &&
         categoryVariable_->group() == "MetaData";
}

}  // namespace ufo
