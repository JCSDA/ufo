/* -----------------------------------------------------------------------------
 * (C) British Crown Copyright 2026 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>

#include "ufo/filters/RecordThresholdRejection.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"

#include "oops/util/Logger.h"

#include "ufo/filters/getScalarOrFilterData.h"

namespace ufo {

constexpr char RecordThresholdRejectionTypeParameterTraitsHelper::enumTypeName[];
constexpr util::NamedEnumerator<RecordThresholdRejectionType>
  RecordThresholdRejectionTypeParameterTraitsHelper::namedValues[];

constexpr char RecordThresholdRejectionDataOrderParameterTraitsHelper::enumTypeName[];
constexpr util::NamedEnumerator<RecordThresholdRejectionDataOrder>
  RecordThresholdRejectionDataOrderParameterTraitsHelper::namedValues[];

// -----------------------------------------------------------------------------

RecordThresholdRejection::RecordThresholdRejection(ioda::ObsSpace & obsdb,
                                   const Parameters_ & parameters,
                                   ioda::ObsDataVector<int> & flags,
                                   ioda::ObsDataVector<float> & obserr)
  : FilterBase(obsdb, parameters, flags, obserr), parameters_(parameters)
{
  oops::Log::trace() << "RecordThresholdRejection constructor" << std::endl;
}

// -----------------------------------------------------------------------------

RecordThresholdRejection::~RecordThresholdRejection() {
  oops::Log::trace() << "RecordThresholdRejection destructor" << std::endl;
}

// -----------------------------------------------------------------------------

/// \brief Filter to perform actions for all entries within a record after a
/*! certain value is reached (in ascending or descending data ordering). 
 * See RecordThresholdRejectionParameters for the documentation of the parameters 
 * controlling this filter.
 */
void RecordThresholdRejection::applyFilter(const std::vector<bool> & apply,
                                   const Variables & filtervars,
                                   std::vector<std::vector<bool>> & flagged) const {
  oops::Log::trace() << "RecordThresholdRejection applyFilter start" << std::endl;

  if (obsdb_.obs_group_vars().empty()) {
    std::ostringstream msg;
    msg << "Group variables configuration is empty, variables are not grouped into records";
    throw eckit::UserError(msg.str(), Here());
  }

  std::vector<float> thrValue(obsdb_.nlocs());
  thrValue = getScalarOrFilterData(parameters_.threshold_value.value(), data_);

  std::vector<float> thrVariable(obsdb_.nlocs());
  thrVariable = getScalarOrFilterData(parameters_.threshold_variable.value(), data_);

  RecordThresholdRejectionDataOrder dataOrder = parameters_.data_order;

  // True channel count from ObsSpace (0 for non-radiance data).
  const size_t nChans = obsdb_.nchans();
  // Effective channel count for indexing: use one slot for non-channel data.
  const size_t nEffectiveChans = std::max(nChans, size_t{1});

  size_t nActualVars = filtervars.nvars();
  if (nChans > 0) {
    // For multi-channel data "nvars" is the number of simulated variables times
    // the number of channels.
    if (filtervars.nvars() % nChans != 0)
      throw eckit::BadValue("filtervars.nvars() is not divisible by obsdb_.nchans()", Here());
    nActualVars = filtervars.nvars() / nChans;
  }

  const float missing = util::missingValue<float>();

  const std::vector<size_t> & recordNumbers = obsdb_.recidx_all_recnums();

  typedef bool (*Predicate)(float, float);
  Predicate inequality = nullptr;
  switch (parameters_.rejection_type) {
    case RecordThresholdRejectionType::LT:
      inequality = [](float a, float b) { return a < b; };
      break;
    case RecordThresholdRejectionType::LTE:
      inequality = [](float a, float b) { return a <= b; };
      break;
    case RecordThresholdRejectionType::GT:
      inequality = [](float a, float b) { return a > b; };
      break;
    case RecordThresholdRejectionType::GTE:
      inequality = [](float a, float b) { return a >= b; };
      break;
    default:
      // should never get here because the parameter is validated by the ParameterTraits,
      // but just in case
      throw eckit::UserError("Unknown rejection type", Here());
  }

  // Loop over the number of actual variables
  for (size_t iVar=0; iVar < nActualVars; ++iVar) {
    // Loop over the unique records
    for (size_t iRecord : recordNumbers) {
      std::vector<size_t> recordIdxs = obsdb_.recidx_vector(iRecord);
      if (dataOrder == RecordThresholdRejectionDataOrder::Descending) {
        std::reverse(recordIdxs.begin(), recordIdxs.end());
      }

      // For each channel and vertical level count the number of valid observations
      for (size_t iChan=0; iChan < nEffectiveChans; ++iChan) {
        const size_t iFilterVar = iVar * nEffectiveChans + iChan;
        bool reject = false;

        // Count the number of valid observations in this record
        for (size_t jobs : recordIdxs) {
          if (apply[jobs]) {
              if  (thrVariable[jobs] != missing && thrValue[jobs] != missing
                   && inequality(thrVariable[jobs], thrValue[jobs])) {
                reject = true;
              }
            if (reject) flagged[iFilterVar][jobs] = true;
          }
        }
      }
    }
  }
  oops::Log::trace() << "RecordThresholdRejection applyFilter complete" << std::endl;
}

// -----------------------------------------------------------------------------

void RecordThresholdRejection::print(std::ostream & os) const {
  os << "RecordThresholdRejection: config = " << parameters_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
