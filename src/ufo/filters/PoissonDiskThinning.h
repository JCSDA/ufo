/*
 * (C) Copyright 2019 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_POISSONDISKTHINNING_H_
#define UFO_FILTERS_POISSONDISKTHINNING_H_

#include <array>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include <boost/optional.hpp>

#include "ioda/ObsDataVector.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/filters/FilterBase.h"
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
  template <typename Key, typename Value> class ScalarOrMap;
}

namespace ufo {

class DistanceCalculator;
class EquispacedBinSelector;
class ObsAccessor;
class PoissonDiskThinningParameters;
class RecursiveSplitter;
class SpatialBinSelector;

/// \brief Thins observations by iterating over them in random order and retaining each observation
/// lying outside the _exclusion volumes_ (ellipsoids or cylinders) surrounding observations that
/// have already been retained.
///
/// See PoissonDiskThinningParameters for the documentation of the available options.
class PoissonDiskThinning : public FilterBase,
                          private util::ObjectCounter<PoissonDiskThinning> {
 public:
  static const std::string classname() {return "ufo::PoissonDiskThinning";}

  PoissonDiskThinning(ioda::ObsSpace &obsdb, const eckit::Configuration &config,
                      std::shared_ptr<ioda::ObsDataVector<int> > flags,
                      std::shared_ptr<ioda::ObsDataVector<float> > obserr);

  ~PoissonDiskThinning() override;

 private:
  struct ObsData;

  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {return QCflags::thinned;}

  ObsAccessor createObsAccessor() const;

  /// \brief Collect all observation data components used for thinning.
  ///
  /// \param[out] numSpatialDims
  ///   Number of spatial dimensions used for thinning (3 if thinning by latitude and longitude,
  ///   0 otherwise).
  /// \param[out] numNonspatialDims
  ///   Number of non-spatial dimensions used for thinning.
  ObsData getObsData(const ObsAccessor &obsAccessor,
                     int &numSpatialDims, int &numNonspatialDims) const;

  template <typename ValueType>
  void validateSpacings(const util::ScalarOrMap<int, ValueType> &spacingsByPriority,
      const std::string &parameterName) const;

  std::vector<size_t> getValidObservationIds(const std::vector<bool> &apply,
                                             const ObsAccessor &obsAccessor) const;

  /// Initialise random number generators used for shuffling with the same seed on
  /// all processes belonging to the given communicator.
  void synchroniseRandomNumberGenerators(const eckit::mpi::Comm &comm) const;

  void groupObservationsByPriority(const std::vector<size_t> &validObsIds,
                                   const ObsAccessor &obsAccessor,
                                   RecursiveSplitter &splitter) const;

  void thinCategory(const ObsData &obsData,
                    const std::vector<size_t> &obsIdsInCategory,
                    const RecursiveSplitter &prioritySplitter,
                    int numSpatialDims,
                    int numNonspatialDims,
                    std::vector<bool> &isThinned) const;

  /// Thin observations belonging to a single category.
  template <int numDims>
  void thinCategory(const ObsData &obsData,
                    const std::vector<size_t> &obsIdsInCategory,
                    const RecursiveSplitter &prioritySplitter,
                    int numSpatialDims,
                    std::vector<bool> &isThinned) const;

  template <int numDims>
  std::array<float, numDims> getObservationPosition(
      size_t obsId, const ObsData &obsData) const;

  template <int numDims>
  std::array<float, numDims> getExclusionVolumeSemiAxes(
      size_t obsId, const ObsData &obsData) const;

 private:
  std::unique_ptr<PoissonDiskThinningParameters> options_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_POISSONDISKTHINNING_H_
