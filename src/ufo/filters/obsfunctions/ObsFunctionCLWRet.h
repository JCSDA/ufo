/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONCLWRET_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONCLWRET_H_

#include <string>
#include <vector>

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

///
/// \brief Retrieve cloud liquid water from AMSU-A 23.8 GHz and 31.4 GHz channels.
///
/// Reference: Grody et al. (2001)
///
/// Determination of precipitable water and cloud liquid water over oceans from
/// the NOAA 15 advanced microwave sounding unit.
///

class ObsFunctionCLWRet : public ObsFunctionBase {
 public:
  explicit ObsFunctionCLWRet(const eckit::LocalConfiguration conf
                                       = eckit::LocalConfiguration());
  ~ObsFunctionCLWRet();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
  inline static float getBadValue() {return bad_clwret_value_;}
 private:
  static constexpr float bad_clwret_value_ = 1000.f;
  ufo::Variables invars_;
  std::vector<std::string> group_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONCLWRET_H_
