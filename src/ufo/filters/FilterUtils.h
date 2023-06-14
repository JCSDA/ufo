/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_FILTERUTILS_H_
#define UFO_FILTERS_FILTERUTILS_H_

#include <memory>
#include <string>
#include <vector>

namespace ioda {
template <typename DATATYPE> class ObsDataVector;
}

namespace ufo {

class Variables;

enum class UnselectLocationIf {
  ANY_FILTER_VARIABLE_REJECTED,
  ALL_FILTER_VARIABLES_REJECTED
};

/// \brief Unselect locations at which observations of any filter variable or all filter variables
/// (depending on the `mode` parameter) have been rejected.
///
/// \param[inout] selected
///   A Boolean vector with as many elements as there are locations held by the calling process
///   (typically the `apply` vector produced by a filter's `where` clause).
///   On output, elements of this vector are set to false at all locations where any (if `mode` is
///   `ANY_FILTER_VARIABLE_REJECTED`) or all (if `mode` is `ALL_FILTER_VARIABLE_REJECTED`) of the
///   variables `filtervars` have been rejected.
/// \param[in] filtervars
///   A subset of simulated variables.
/// \param[in] qcflags
///   QC flags of all simulated variables.
/// \param[in] mode
///   Indicates whether to unselect locations where any or all of the variables `filtervars` have
///   been rejected.
void unselectRejectedLocations(std::vector<bool> &selected,
                               const ufo::Variables &filtervars,
                               const ioda::ObsDataVector<int> &qcflags,
                               UnselectLocationIf mode,
                               const std::vector<size_t> &obs_inds = {});

}  // namespace ufo

#endif  // UFO_FILTERS_FILTERUTILS_H_
