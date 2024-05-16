/*
 * (C) Crown copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_RECORDNUMBERTOOBSSPACE_H_
#define UFO_FILTERS_OBSFUNCTIONS_RECORDNUMBERTOOBSSPACE_H_

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variable.h"
#include "ufo/filters/Variables.h"

namespace ufo {

/// \brief Write record number to a variable which can be saved to the ObsSpace.
/// This ObsFunction can only be used for data that have been grouped into records.
class RecordNumberToObsSpace : public ObsFunctionBase<int> {
 public:
  explicit RecordNumberToObsSpace(const eckit::LocalConfiguration &);
  ~RecordNumberToObsSpace();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<int> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_RECORDNUMBERTOOBSSPACE_H_
