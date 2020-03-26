/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONOBSERRORMEAN_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONOBSERRORMEAN_H_

#include <string>
#include <vector>

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

///
/// \brief Calculate symmetric (mean) observation error from the average of the observed and
/// simulated cloud amount.
///
/// References:
/// (1) Geer, A. J. and P. Bauer (2011). Observation errors in all-sky data
///     assimilation. Quart. J. Roy. Meteorol. Soc. 137, 2024–2037.
/// (2) Geer, A. J., P. Bauer, and S. J. English (2012). Assimilating AMSU-A
///     temperature sounding channels in the presence of cloud and precipitation.
///     Published simultaneously as ECMWF Technical Memoranda 670 and ECMWF/EUMETSAT
///     fellowship reports 24.
///

class ObsFunctionObsErrorMean : public ObsFunctionBase {
 public:
  explicit ObsFunctionObsErrorMean(const eckit::LocalConfiguration);
  ~ObsFunctionObsErrorMean();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  std::vector<int> channels_;
  const eckit::LocalConfiguration conf_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONOBSERRORMEAN_H_
