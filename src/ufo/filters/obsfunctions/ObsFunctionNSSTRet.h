/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONNSSTRET_H_
#define UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONNSSTRET_H_

#include <string>
#include <vector>

#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

// -----------------------------------------------------------------------------
// QC using retrieved near sea surface temperatur (NSST) from radiances
// 2-step QC procedures:
// (1) Perform NSST retrieval from radiances at obs location, and obtained
//     increment of NSST from its first guess value
// (2) For surface sensitive channels, remove them from assimilation if the
//     retrieved NSST increment from step (1) is larger than a pre-defined
//     threshold
// Output:
// 0 = channel is retained for assimilation
// 1 = channel is not retained for assimilation
// -----------------------------------------------------------------------------

class ObsFunctionNSSTRet : public ObsFunctionBase {
 public:
  explicit ObsFunctionNSSTRet(const eckit::LocalConfiguration);
  ~ObsFunctionNSSTRet();

  void compute(const ObsFilterData &,
               ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  std::string group_;
  std::vector<int> channels_;
  const eckit::LocalConfiguration conf_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_OBSFUNCTIONNSSTRET_H_
