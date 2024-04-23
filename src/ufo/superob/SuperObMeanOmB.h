/*
 * (C) Crown copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_SUPEROB_SUPEROBMEANOMB_H_
#define UFO_SUPEROB_SUPEROBMEANOMB_H_

#include <vector>

#include "ufo/superob/SuperObBase.h"

namespace ufo {

/// \brief Simple superob computation.
///
/// For each record consisting of N observations O and background values B,
/// a superob S is computed as follows:
///
/// S = B_l + sum_{k = 1}^{N} (O_k - B_k) / N
///
/// where B_l is the background value at the superob location.
/// In this algorithm the superob location is the first location in the record
/// with a passing QC flag.
class SuperObMeanOmB : public SuperObBase {
 public:
  explicit SuperObMeanOmB(const GenericSuperObParameters &,
                          const ObsFilterData &,
                          const std::vector<bool> &,
                          const Variables &,
                          const ioda::ObsDataVector<int> &,
                          std::vector<std::vector<bool>> &);
  ~SuperObMeanOmB() {}

 private:
  /// Compute superobs for each record.
  void computeSuperOb(const std::vector<std::size_t> &,
                      const std::vector<float> &,
                      const std::vector<float> &,
                      const ioda::ObsDataRow<int> &,
                      std::vector<float> &,
                      std::vector<bool> &) const override;
};

}  // namespace ufo

#endif  // UFO_SUPEROB_SUPEROBMEANOMB_H_
