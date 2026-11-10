/*
 * (C) Copyright 2026 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/PropagateOneValuePerRecord.h"

#include <cstddef>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/GeoVaLs.h"

namespace ufo {

static ObsFunctionMaker<PropagateOneValuePerRecord<float>>
                       makerPropagateOneValuePerRecordFloat_("PropagateOneValuePerRecord");
static ObsFunctionMaker<PropagateOneValuePerRecord<int>>
                       makerPropagateOneValuePerRecordInt_("PropagateOneValuePerRecord");
static ObsFunctionMaker<PropagateOneValuePerRecord<std::string>>
                       makerPropagateOneValuePerRecordString_("PropagateOneValuePerRecord");
static ObsFunctionMaker<PropagateOneValuePerRecord<util::DateTime>>
                       makerPropagateOneValuePerRecordDateTime_("PropagateOneValuePerRecord");
// -----------------------------------------------------------------------------

template <typename PropagateVariableType>
PropagateOneValuePerRecord<PropagateVariableType>::PropagateOneValuePerRecord(
        const eckit::LocalConfiguration & conf)
  : invars_() {
  oops::Log::trace() << "PropagateOneValuePerRecord constructor" << std::endl;
  // Validate and deserialize options
  options_.validateAndDeserialize(conf);

  invars_ += options_.variableToPropagate.value();

  if (options_.propagateToCompanionAsWell.value() == true) {
    if (!options_.propagationDirection.value()) {
      throw eckit::UserError(
         "Propagation direction must be specified when propagating to companion",
         Here());
    }
  } else {
    if (options_.propagationDirection.value()) {
      throw eckit::UserError(
         "Propagation direction must not be specified if not propagating to companion",
         Here());
    }
  }
}

// -----------------------------------------------------------------------------
template <typename PropagateVariableType>
PropagateOneValuePerRecord<PropagateVariableType>::~PropagateOneValuePerRecord() {
  oops::Log::trace() << "PropagateOneValuePerRecord destructor" << std::endl;
}

// -----------------------------------------------------------------------------
template <typename PropagateVariableType>
void PropagateOneValuePerRecord<PropagateVariableType>::compute(const ObsFilterData & in,
                                ioda::ObsDataVector<PropagateVariableType> & out) const {
  oops::Log::trace() << "PropagateOneValuePerRecord compute start" << std::endl;

  // ObsSpace.
  ioda::ObsSpace & obsdb = in.obsspace();

  // Ensure observations have been grouped into records.
  if (obsdb.obs_group_vars().empty())
    throw eckit::UserError("Group variables configuration is empty", Here());

  // Correspondence between record numbers and indices in the data sample.
  const std::vector<std::size_t> &recnums = obsdb.recidx_all_recnums();

  std::size_t nrecsToPropagate;    // number of records to iterate over
  if (options_.propagateToCompanionAsWell.value() == false) {
    nrecsToPropagate = recnums.size();

  } else {     // check for companion records
    // Check the ObsSpace has been extended. If this is not the case
    // then it will not be possible to access records in the original and
    // extended sections of the ObsSpace.
    if (!obsdb.has("MetaData", "extendedObsSpace"))
      throw eckit::UserError("The extended obs space has not been produced", Here());
    // Number of records in the original ObsSpace (extending obs space always
    // doubles number of records).
    nrecsToPropagate = recnums.size() / 2;
  }

  const PropagateVariableType missing = util::missingValue<PropagateVariableType>();

  ioda::ObsDataVector<PropagateVariableType> varin(obsdb,
      options_.variableToPropagate.value().toOopsObsVariables());

  in.get(options_.variableToPropagate.value(), varin, options_.skipDerived.value());
  if (varin.varnames().channels().size() != out.varnames().channels().size()) {
    throw eckit::UserError("Number of channels in the variable to propagate does not match "
                           "the number of channels in the output variable", Here());
  }

  for (std::size_t jrec = 0; jrec < nrecsToPropagate; ++jrec) {
    std::size_t jToPropagateFrom;  // index for the record from which to propagate from
    if (options_.propagateToCompanionAsWell.value() ==
        false ||
        options_.propagationDirection.value() ==
        PropagateOneValuePerRecordPropagationDirection::OriginalToExtended) {
      jToPropagateFrom = jrec;

    } else {    // PropagateOneValuePerRecordPropagationDirection::ExtendedToOriginal
      // Corresponding extended obs space record is `recnums.size() / 2` away
      jToPropagateFrom = jrec + (recnums.size() / 2);
    }

    const std::vector<std::size_t> &locs =
      obsdb.recidx_vector(recnums[jToPropagateFrom]);

    std::optional<std::size_t> locIndex
      = findLocIndex(options_.inRecordLocationOfValueToPropagate.value(),
                     locs.size());

    for (size_t iChan = 0; iChan < out.nvars(); ++iChan) {
      PropagateVariableType valueToPropagate;
      if (locs.empty() || !locIndex.has_value()) {
        valueToPropagate = missing;
      } else {
        valueToPropagate = varin[iChan][locs[locIndex.value()]];
      }

      for (std::size_t loc : locs)
        out[iChan][loc] = valueToPropagate;

      if (options_.propagateToCompanionAsWell.value() == true) {
        std::size_t companionRecordIndex
          = (jToPropagateFrom + (recnums.size() / 2)) % recnums.size();

        const std::vector<std::size_t> &locsCompanion =
          obsdb.recidx_vector(recnums[companionRecordIndex]);

        for (std::size_t loc : locsCompanion)
          out[iChan][loc] = valueToPropagate;
      }
    }
  }
  oops::Log::trace() << "PropagateOneValuePerRecord compute complete" << std::endl;
}

// -----------------------------------------------------------------------------
template <typename PropagateVariableType>
std::optional<std::size_t> PropagateOneValuePerRecord<PropagateVariableType>::findLocIndex(
                                    const int value, const std::size_t locsSize) const {
  int output;
  if (value < 0) {
    output = static_cast<int>(locsSize) + value;
  } else {
    output = value;
  }
  if (output < 0 || output >= static_cast<int>(locsSize)) {
    return std::nullopt;
  }
  return static_cast<std::size_t>(output);
}

// -----------------------------------------------------------------------------
template <typename PropagateVariableType>
const ufo::Variables &
    PropagateOneValuePerRecord<PropagateVariableType>::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
