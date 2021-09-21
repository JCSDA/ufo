/*
 * (C) Copyright 2017-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_LOCATIONS_H_
#define UFO_LOCATIONS_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "eckit/mpi/Comm.h"
#include "ioda/distribution/Distribution.h"
#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace ufo {

/// \brief Locations class to handle simple lat-lon-time locations
class Locations : public util::Printable,
                  private util::ObjectCounter<Locations> {
 public:
  static const std::string classname() {return "ufo::Locations";}

  /// \brief constructor from passed \p lons, \p lats, \p times
  Locations(const std::vector<float> & lons, const std::vector<float> & lats,
            const std::vector<util::DateTime> & times,
            std::shared_ptr<const ioda::Distribution>);
  /// \brief constructor used in oops tests
  Locations(const eckit::Configuration &, const eckit::mpi::Comm &);

  /// append locations with more locations
  Locations & operator+=(const Locations &);

  /// find which observations are in the (\p t1, \p t2] time window
  std::vector<bool> isInTimeWindow(const util::DateTime & t1, const util::DateTime & t2) const;

  /// size of locations
  size_t size() const;

  /// accessor to the observations MPI distribution
  const std::shared_ptr<const ioda::Distribution> & distribution() const {return dist_;}
  /// accessor to observation longitudes (on current MPI task)
  const std::vector<float> & lons() const {return lons_;}
  /// accessor to observation latitudes (on current MPI task)
  const std::vector<float> & lats() const {return lats_;}
  /// accessor to DateTimes (on current MPI task)
  const std::vector<util::DateTime> & times() const {return times_;}


 private:
  void print(std::ostream & os) const override;

  std::shared_ptr<const ioda::Distribution> dist_;   /// observations MPI distribution
  std::vector<float> lons_;            /// longitudes on current MPI task
  std::vector<float> lats_;            /// latitudes on current MPI task
  std::vector<util::DateTime> times_;  /// times of observations on current MPI task
};

}  // namespace ufo

#endif  // UFO_LOCATIONS_H_
