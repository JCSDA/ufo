/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_SAMPLEDLOCATIONS_H_
#define UFO_SAMPLEDLOCATIONS_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "eckit/mpi/Comm.h"

#include "ioda/distribution/Distribution.h"
#include "ioda/ObsGroup.h"

#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Range.h"

namespace eckit {
  class Configuration;
}

namespace ufo {

/// \brief Collection of interpolation paths sampling the observation locations.
///
/// The location of a single observation may be an extended region of space and be sampled by
/// more than one path. Conversely, some paths may sample the locations of multiple
/// observations.
///
/// \note Currently each interpolation path is assumed to be vertical. There are plans to add
/// support for slanted paths in the future.
class SampledLocations : public util::Printable,
                         private util::ObjectCounter<SampledLocations> {
 public:
  static const std::string classname() {return "ufo::SampledLocations";}

  /// \brief Create a collection of vertical paths with specified longitudes, latitudes and times.
  ///
  /// \param lons
  ///   Interpolation path longitudes.
  /// \param lats
  ///   Interpolation path latitudes.
  /// \param times
  ///   Interpolation times.
  /// \param dist
  ///   MPI distribution of the observation locations sampled by the newly constructed paths.
  /// \param pathsGroupedByLocations
  ///   (Optional) A vector mapping the index of each observation location to the range of
  ///   indices of the paths sampling that location. Specifically, the ith location is deemed to
  ///   be sampled by the paths with indices ranging from `pathsGroupedByLocations[i].begin` up to
  ///   but not including `pathsGroupedByLocations[i].end`.
  ///
  ///   This vector may be left empty if each location is sampled solely by the path with the same
  ///   index.
  SampledLocations(const std::vector<float> & lons, const std::vector<float> & lats,
                   const std::vector<util::DateTime> & times,
                   std::shared_ptr<const ioda::Distribution> dist,
                   std::vector<util::Range<size_t>> pathsGroupedByLocation = {});

  /// \brief constructor used in oops tests
  SampledLocations(const eckit::Configuration &, const eckit::mpi::Comm &);

  /// append paths from another collection
  ///
  /// TODO(wsmigaj): Check if this method can be removed, especially since it doesn't update dist_
  /// even though it should.
  SampledLocations & operator+=(const SampledLocations &);

  /// find which observations are in the (\p t1, \p t2] time window
  std::vector<bool> isInTimeWindow(const util::DateTime & t1, const util::DateTime & t2) const;

  /// number of paths belonging to this collection
  size_t size() const;

  /// accessor to the observation locations' MPI distribution
  const std::shared_ptr<const ioda::Distribution> & distribution() const {return dist_;}
  /// accessor to the longitudes of successive interpolation paths (on current MPI task)
  /// TO BE CONSOLIDATED
  std::vector<float> lons() const;
  const std::vector<double> & longitudes() const {return lons_;}
  /// accessor to the latitudes of successive interpolation paths (on current MPI task)
  /// TO BE CONSOLIDATED
  std::vector<float> lats() const;
  const std::vector<double> & latitudes() const {return lats_;}
  /// accessor to the DateTimes of successive interpolation paths (on current MPI task)
  const std::vector<util::DateTime> & times() const {return times_;}

  /// \brief Returns a vector mapping the index of each observation location to the range of
  /// indices of the paths sampling that location.
  ///
  /// If the returned vector `v` is empty, it means each location is sampled solely by the path
  /// with the same index. Otherwise the `i`th location is sampled by the paths with indices
  /// from `v[i].begin` up to but not including `v[i].end`.
  const std::vector<util::Range<size_t>> &pathsGroupedByLocation() const {
    return pathsGroupedByLocation_;
  }

  /// \brief Number of locations sampled by this collection of paths.
  size_t nlocs() const;

  /// \brief Returns true if and only if all locations are sampled exactly once, in ascending order.
  ///
  /// Another way to put it is that each location is sampled solely by the interpolation path with
  /// the same index.
  bool areLocationsSampledOnceAndInOrder() const;

 private:
  void initializeObsGroup(size_t npaths);
  void print(std::ostream & os) const override;

  std::shared_ptr<const ioda::Distribution> dist_;   /// sampled locations' MPI distribution
  ioda::ObsGroup og_;  /// interpolation paths on current MPI task
  // Note on implementation: at time of writing, the ObsGroup does not support the util::DateTime
  // type, so instead of doing expensive conversions of times to/from string representations, we
  // opt to keep the times outside the ObsGroup.
  std::vector<util::DateTime> times_;  /// interpolation path times on current MPI task
  std::vector<double> lons_;  /// interpolation path longitudes on current MPI task
  std::vector<double> lats_;  /// interpolation path latitudes on current MPI task
  std::vector<util::Range<size_t>> pathsGroupedByLocation_;
};

}  // namespace ufo

#endif  // UFO_SAMPLEDLOCATIONS_H_
