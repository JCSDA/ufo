/*
 * (C) Copyright 2022 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsFunctionSelectStatistic.h"

#include <algorithm>
#include <numeric>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/utils/RecursiveSplitter.h"

namespace ufo {

static ObsFunctionMaker<SelectStatistic> makerSelectStatistic_("SelectStatistic");

// -----------------------------------------------------------------------------

SelectStatistic::SelectStatistic(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Check options
  options_.validateAndDeserialize(conf);

  // Create variable and add to invars_
  for (const Variable & var : options_.variable.value()) {
    invars_ += var;
  }
  ASSERT(invars_.size() == 1);  // only works with one (possibly multi-channel) x variable
}

// -----------------------------------------------------------------------------

SelectStatistic::~SelectStatistic() {}

// -----------------------------------------------------------------------------

ufo::ObsAccessor SelectStatistic::createObsAccessor(const ioda::ObsSpace &obsdb) const {
  if (obsdb.obs_group_vars().empty()) {
    return ObsAccessor::toAllObservations(obsdb);
  } else {
    return ObsAccessor::toObservationsSplitIntoIndependentGroupsByRecordId(obsdb);
  }
}

// -----------------------------------------------------------------------------

// template <typename InputType>
void SelectStatistic::compute(const ObsFilterData & in, ioda::ObsDataVector<int> & out) const {
  // initialize
  out.zero();
  size_t nchans = out.nvars();

  // ObsSpace.
  ioda::ObsSpace & obsdb = in.obsspace();

  // Vector of locations that pass the 'where' clause in the sample
  // (all true if there is no where clause).
  const std::vector<bool> apply = processWhere(options_.where, in, options_.whereOperator);

  const float missing = util::missingValue<float>();
  ObsAccessor obsAccessor = createObsAccessor(obsdb);
  const std::vector<size_t> validObsIds = obsAccessor.getValidObservationIds(apply);
  oops::Log::debug() << "validObsIds: " << validObsIds << std::endl;
  RecursiveSplitter splitter = obsAccessor.splitObservationsIntoIndependentGroups(validObsIds);
  ioda::ObsDataVector<float> varin(in.obsspace(), invars_[0].toOopsObsVariables());
  in.get(invars_[0], varin);
  oops::Log::debug() << "out.nvars(): " << out.nvars() << std::endl;

  for (ufo::RecursiveSplitter::Group recordGroup : splitter.groups()) {
    // Get input vector of valid values in channel ichan:
    for (size_t ichan = 0; ichan < nchans; ++ichan) {
      oops::Log::debug() << "ichan: " << ichan << std::endl;
      std::vector<float> inputVec;
      std::vector<float> obs_indices;
      // std::vector<float> all_inputVec;
      std::vector<float> all_obs_indices;
      // Loop over obs within a profile:
      for (size_t validObsIndex : recordGroup) {
        all_obs_indices.push_back(validObsIndex);  // include missing vals
        if (varin[ichan][validObsIds[validObsIndex]] != missing) {
          inputVec.push_back(varin[ichan][validObsIds[validObsIndex]]);
          obs_indices.push_back(validObsIndex);
        }
      }
      if (inputVec.size() > 0) {  // leave output as 0's if input var all missing for whole record
        if (options_.selectMax) {
          oops::Log::debug() << "Select max." << std::endl;
          const size_t maxInd = std::distance(inputVec.begin(),
                                              std::max_element(inputVec.begin(), inputVec.end()));
          out[ichan][validObsIds[obs_indices[maxInd]]] = 1;
        }
        if (options_.selectMin) {
          oops::Log::debug() << "Select min." << std::endl;
          const size_t minInd = std::distance(inputVec.begin(),
                                              std::min_element(inputVec.begin(), inputVec.end()));
          out[ichan][validObsIds[obs_indices[minInd]]] = 1;
        }
        if (options_.selectMean) {
          oops::Log::debug() << "Select mean." << std::endl;
          size_t inputSize = inputVec.size();
          const float inputMean = std::accumulate(inputVec.begin(), inputVec.end(), 0.0) /
                                                  inputSize;
          auto i = std::min_element(inputVec.begin(), inputVec.end(), [=] (float x, float y)
          {
              return abs(x - inputMean) < abs(y - inputMean);
          });
          const size_t meanInd = std::distance(inputVec.begin(), i);
          out[ichan][validObsIds[obs_indices[meanInd]]] = 1;
        }
        if (options_.selectMedian) {
          oops::Log::debug() << "Select median." << std::endl;
          size_t inputSize = inputVec.size();
          std::vector<size_t> idx(inputSize);
          std::iota(idx.begin(), idx.end(), 0);
          std::stable_sort(idx.begin(), idx.end(), [&inputVec](size_t i1, size_t i2)
                                                   {return inputVec[i1] < inputVec[i2];});
          const size_t medianInd = idx[(inputSize-1)/2];
          out[ichan][validObsIds[obs_indices[medianInd]]] = 1;
        }
      } else {  // inputVec all missing
        if (options_.forceSelect) {
          oops::Log::debug() <<
              "All missing values in this record, but force to select one obs." << std::endl;
          out[ichan][validObsIds[all_obs_indices[0]]] = 1;
        }
      }  // record all missing values or not
    }  // for each channel
  }  // for each record
}

// -----------------------------------------------------------------------------

const ufo::Variables & SelectStatistic::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
