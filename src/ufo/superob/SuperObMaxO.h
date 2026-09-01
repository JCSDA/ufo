/*
 * (C) Crown copyright 2025, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_SUPEROB_SUPEROBMAXO_H_
#define UFO_SUPEROB_SUPEROBMAXO_H_

#include <limits>
#include <vector>

#include "ufo/superob/SuperObBase.h"

namespace ufo {

/// \brief Compute the maximum non-missing observation in a record.
///
/// For each record consisting of N observations O, a superob S is computed as
/// follows:
///
/// S = max_{k = 1}^{N}(O_k), where O_k != missing.
///
/// Furthermore, O_k which are not flagged as passing QC or as passive are
/// excluded from the computation.
///
/// In this algorithm the superob location is the first location in the record
/// with a passing QC flag.
class SuperObMaxO : public SuperObBase {
 public:
  using Parameters_ = GenericSuperObParameters;

  explicit SuperObMaxO(const Parameters_ &, const ObsFilterData &,
                       const std::vector<bool> &, const Variables &,
                       const ioda::ObsDataVector<int> &,
                       std::vector<std::vector<bool>> &);
  ~SuperObMaxO() {}

  bool requireHofX() const override { return false; }

 private:
  /// Compute superobs for each record. Returns `true` if a
  /// superob was successfully computed for the record, `false` otherwise.
  bool computeSuperOb(const std::vector<std::size_t> &,
                      const std::vector<float> &, const std::vector<float> &,
                      const ioda::ObsDataRow<int> &, std::vector<float> &,
                      std::vector<bool> &) const override;

  const Parameters_ &params_;
};

}  // namespace ufo

#endif  // UFO_SUPEROB_SUPEROBMAXO_H_
