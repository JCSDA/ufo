/*
 * (C) Crown copyright 2025, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_USENEARESTNEIGHBORS_H_
#define UFO_FILTERS_USENEARESTNEIGHBORS_H_

#include <memory>
#include <string>
#include <vector>

#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"
#include "ufo/filters/UseNearestNeighborsParameters.h"
#include "ufo/usenearestneighbors/UseNearestNeighborsAlgorithmBase.h"

namespace ufo {

/// \brief Filter that uses nearest-neighbor information to perform operations
/// on observations relative to their spatial neighbors.
///
/// The filter is designed to work with an obs space that has already been
/// populated (typically by the Find Nearest Neighbors filter) with variables
/// identifying, for each query-point location, the first, second, third, etc.
/// nearest reference observation locations. The specific operation is
/// controlled by the \c algorithm block; currently supported algorithms are:
///
/// - "gather and match timestamp": for each query point, looks up a value
///   from the nearest reference observations using an exact match on identifier
///   and timestamp.
/// - "reference point variables mean": for each query point, averages a
///   value gathered from the nearest reference observations, grouped by an
///   integer index value.
/// - "local plane fit": interpolates/extrapolates a scalar from the nearest
///   reference observations using weighted least-squares plane fitting, with
///   fallback to inverse distance weighting.
class UseNearestNeighbors : public FilterBase,
                            private util::ObjectCounter<UseNearestNeighbors> {
 public:
  typedef UseNearestNeighborsParameters Parameters_;
  static const std::string classname() { return "ufo::UseNearestNeighbors"; }

  UseNearestNeighbors(ioda::ObsSpace &, const Parameters_ &,
                      ioda::ObsDataVector<int> &,
                      ioda::ObsDataVector<float> &);
  ~UseNearestNeighbors();

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override { return QCflags::pass; }
  Parameters_ options_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_USENEARESTNEIGHBORS_H_
