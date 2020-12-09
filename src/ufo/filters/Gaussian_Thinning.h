/*
 * (C) Copyright 2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_GAUSSIAN_THINNING_H_
#define UFO_FILTERS_GAUSSIAN_THINNING_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include <boost/optional.hpp>

#include "ioda/ObsDataVector.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/GaussianThinningParameters.h"
#include "ufo/filters/QCflags.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace util {
  class DateTime;
}

namespace ufo {

class DistanceCalculator;
class EquispacedBinSelector;
class GaussianThinningParameters;
class ObsAccessor;
class RecursiveSplitter;
class SpatialBinSelector;

/// \brief Group observations into grid cells and preserve only one observation in each cell.
///
/// Cell assignment can be based on an arbitrary combination of:
/// - horizontal position
/// - vertical position (in terms of air pressure)
/// - time
/// - category (arbitrary integer associated with each observation).
///
/// Selection of the observation to preserve in each cell is based on
/// - its position in the cell
/// - optionally, its priority.
///
/// See GaussianThinningParameters for the documentation of the available options.
class Gaussian_Thinning : public FilterBase,
                          private util::ObjectCounter<Gaussian_Thinning> {
 public:
  /// The type of parameters accepted by the constructor of this filter.
  /// This typedef is used by the FilterFactory.
  typedef GaussianThinningParameters Parameters_;

  static const std::string classname() {return "ufo::Gaussian_Thinning";}

  Gaussian_Thinning(ioda::ObsSpace &obsdb, const GaussianThinningParameters &params,
                    std::shared_ptr<ioda::ObsDataVector<int> > flags,
                    std::shared_ptr<ioda::ObsDataVector<float> > obserr);

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::thinned;}

  ObsAccessor createObsAccessor() const;

  void groupObservationsBySpatialLocation(const std::vector<size_t> &validObsIds,
                                          const DistanceCalculator &distanceCalculator,
                                          const ObsAccessor &obsAccessor,
                                          RecursiveSplitter &splitter,
                                          std::vector<float> &distancesToBinCenter) const;

  void groupObservationsByPressure(const std::vector<size_t> &validObsIds,
                                   const DistanceCalculator &distanceCalculator,
                                   const ObsAccessor &obsAccessor,
                                   RecursiveSplitter &splitter,
                                   std::vector<float> &distancesToBinCenter) const;

  void groupObservationsByTime(const std::vector<size_t> &validObsIds,
                               const DistanceCalculator &distanceCalculator,
                               const ObsAccessor &obsAccessor,
                               RecursiveSplitter &splitter,
                               std::vector<float> &distancesToBinCenter) const;

  std::vector<bool> identifyThinnedObservations(
      const std::vector<size_t> &validObsIds,
      const ObsAccessor &obsAccessor,
      const RecursiveSplitter &splitter,
      const std::vector<float> &distancesToBinCenter) const;

  std::function<bool(size_t, size_t)> makeObservationComparator(
      const std::vector<size_t> &validObsIds,
      const std::vector<float> &distancesToBinCenter,
      const ObsAccessor &obsAccessor) const;

  static boost::optional<SpatialBinSelector> makeSpatialBinSelector(
      const GaussianThinningParameters &options);

  static boost::optional<EquispacedBinSelector> makePressureBinSelector(
      const GaussianThinningParameters &options);

  static boost::optional<EquispacedBinSelector> makeTimeBinSelector(
      const GaussianThinningParameters &options, util::DateTime &timeOffset);

  static std::unique_ptr<DistanceCalculator> makeDistanceCalculator(
      const GaussianThinningParameters &options);

 private:
  GaussianThinningParameters options_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_GAUSSIAN_THINNING_H_
