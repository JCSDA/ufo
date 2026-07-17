/*
 * (C) Crown Copyright 2026 UK Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_PERCENTILE_H_
#define UFO_FILTERS_PERCENTILE_H_

#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/PercentileParameters.h"
#include "ufo/filters/QCflags.h"
#include "ufo/utils/RecursiveSplitter.h"

namespace ioda {
template <typename DATATYPE>
class ObsDataVector;
class ObsSpace;
}  // namespace ioda

namespace ufo {

/// \brief Flags observations that fall outside specified percentile ranges.
///
/// This filter selects the central percentile range of data within each group
/// (defined by obs grouping or station_id_variable). Observations outside this
/// range are flagged for rejection.
///
/// Features:
/// - Configurable lower and upper percentiles thresholds
/// - Inclusive or exclusive range comparison
///
/// Percentile threshold calculations use linear interpolation between closest
/// datapoints where necessary to match the requested percentile value. This is
/// equivalent to numpy.percentile with default method='linear'.
class Percentile : public FilterBase, private util::ObjectCounter<Percentile> {
 public:
  typedef PercentileParameters Parameters_;

  static const std::string classname() { return "ufo::Percentile"; }

  Percentile(ioda::ObsSpace &obsdb, const Parameters_ &parameters,
             ioda::ObsDataVector<int> &flags,
             ioda::ObsDataVector<float> &obserr);

  ~Percentile() override;

 private:
  Parameters_ options_;

  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override { return QCflags::percentile; }

  /// Process a single record (station/group) for a single variable
  void processRecord(const RecursiveSplitter::Group &record,
                     const std::vector<size_t> &validObsIds,
                     const std::vector<float> &variableData,
                     const float lowerPercentile, const float upperPercentile,
                     const bool inclusive, std::vector<bool> &isRejected,
                     std::vector<float> &outputData) const;

  /// Calculate percentile value threshold from sorted data, using linear
  /// interpolation between closest datapoints where necessary to match the
  /// requested percentile value. This is equivalent to numpy.percentile with
  /// default method='linear'.
  float calculatePercentileValue(const std::vector<float> &sortedData,
                                 const float percentile) const;
};

}  // namespace ufo

#endif  // UFO_FILTERS_PERCENTILE_H_
