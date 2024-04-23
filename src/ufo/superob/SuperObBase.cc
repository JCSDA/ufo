/*
 * (C) Crown copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/superob/SuperObBase.h"

#include <map>
#include <string>
#include <unordered_map>

#include "oops/util/Logger.h"

namespace ufo {

SuperObBase::SuperObBase(const SuperObParametersBase & params,
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

void SuperObBase::runAlgorithm() const {
  oops::Log::trace() << "SuperObBase::runAlgorithm starting" << std::endl;
  const float missing = util::missingValue<float>();
  const std::size_t nlocs = obsdb_.nlocs();
  const std::vector<std::size_t> & recnums = obsdb_.recidx_all_recnums();

  // Produce lists of locations for each record.
  std::unordered_map<int, std::vector<std::size_t>> locsToUse;

  // Loop over records and fill valid locations in each.
  for (const auto & recnum : recnums) {
    locsToUse[recnum] = {};
    const std::vector<std::size_t> allLocsInRec = obsdb_.recidx_vector(recnum);
    for (size_t jloc : allLocsInRec) {
      if (apply_[jloc]) {
        locsToUse[recnum].push_back(jloc);
      }
    }
  }

  Variables varhofx(filtervars_, "HofX");

  // Loop over each filter variable and compute superobs for each one.
  // Also compute the associated superob errors, and potentially save any
  // auxiliary variables in the algorithm.
  for (size_t jvar = 0; jvar < filtervars_.nvars(); ++jvar) {
    const std::string variableName = filtervars_[jvar].variable();
    std::vector<float> obs(nlocs);
    obsdb_.get_db("ObsValue", variableName, obs);
    std::vector<float> hofx(nlocs);
    data_.get(varhofx.variable(jvar), hofx);
    std::vector<float> superobs(nlocs, missing);
    // Set all entries in `flagged_` to true. One location in each record
    // will be set to `false` in order to indicate where the superob
    // is located.
    std::fill(flagged_[jvar].begin(), flagged_[jvar].end(), true);
    for (const auto & recnum : recnums) {
      computeSuperOb(locsToUse[recnum],
                     obs, hofx, flags_[jvar], superobs, flagged_[jvar]);
    }
    // Save the superob values to the ObsSpace.
    obsdb_.put_db("DerivedObsValue", variableName, superobs);
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

}  // namespace ufo
