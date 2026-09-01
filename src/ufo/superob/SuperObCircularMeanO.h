/*
 * (C) Crown Copyright 2026 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_SUPEROB_SUPEROBCIRCULARMEANO_H_
#define UFO_SUPEROB_SUPEROBCIRCULARMEANO_H_

#include <cmath>
#include <vector>

#include "oops/util/parameters/Parameter.h"
#include "ufo/superob/SuperObBase.h"

namespace ufo {

/// \brief Parameters for the circular mean superobbing algorithm
class SuperObCircularMeanOParameters : public GenericSuperObParameters {
  OOPS_CONCRETE_PARAMETERS(SuperObCircularMeanOParameters,
                           GenericSuperObParameters)
 public:
  /// Lower bound of the range (default: 0)
  oops::Parameter<float> lowerBound{"lower bound", 0.0f, this};

  /// Exclusive upper bound of the range (default: 2π)
  oops::Parameter<float> exclusiveUpperBound{
      "exclusive upper bound", static_cast<float>(2.0 * M_PI), this};
};

/// \brief Compute circular mean observation values in a record.
///
/// For each record consisting of N observations O, a superob S is computed
/// using the circular mean algorithm. This is useful for directional data
/// like wind direction angle. The algorithm follows the scipy implementation:
/// https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.circmean.html
class SuperObCircularMeanO : public SuperObBase {
 public:
  using Parameters_ = SuperObCircularMeanOParameters;

  explicit SuperObCircularMeanO(const SuperObCircularMeanOParameters &,
                                const ObsFilterData &,
                                const std::vector<bool> &, const Variables &,
                                const ioda::ObsDataVector<int> &,
                                std::vector<std::vector<bool>> &);
  ~SuperObCircularMeanO() {}

  bool requireHofX() const override { return false; }

 private:
  /// Compute superobs for each record. Returns `true` if a
  /// superob was successfully computed for the record, `false` otherwise.
  bool computeSuperOb(const std::vector<std::size_t> &,
                      const std::vector<float> &, const std::vector<float> &,
                      const ioda::ObsDataRow<int> &, std::vector<float> &,
                      std::vector<bool> &) const override;

  const SuperObCircularMeanOParameters &params_;
};

}  // namespace ufo

#endif  // UFO_SUPEROB_SUPEROBCIRCULARMEANO_H_
