/*
 * (C) Crown copyright 2021 Met Office. All rights reserved.
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#ifndef UFO_FILTERS_HISTORYCHECK_H_
#define UFO_FILTERS_HISTORYCHECK_H_

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "oops/util/ObjectCounter.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/HistoryCheckParameters.h"
#include "ufo/filters/QCflags.h"
#include "ufo/filters/TrackCheckUtils.h"

namespace eckit {
class Configuration;
}

namespace ioda {
template <typename DATATYPE> class ObsDataVector;
class ObsSpace;
}

namespace ufo {
class ObsAccessor;
class HistoryCheck: public FilterBase,
    private util::ObjectCounter<HistoryCheck> {
 public:
  typedef HistoryCheckParameters Parameters_;

  static const std::string classname() {return "ufo::HistoryCheck";}

///  Track checks and stuck value checks.
///
///  1. Read in wider window of observations
///  2. Apply track check and stuck value check over wider window
///  3. Apply flags over wider window to observations in main window.
  HistoryCheck(ioda::ObsSpace &obsdb, const Parameters_ &parameters,
                 std::shared_ptr<ioda::ObsDataVector<int> > flags,
                 std::shared_ptr<ioda::ObsDataVector<float> > obserr);

  // This constructor is needed for unit testing
  HistoryCheck(ioda::ObsSpace &obsdb, const Parameters_ &parameters,
                 std::shared_ptr<ioda::ObsDataVector<int> > flags,
                 std::shared_ptr<ioda::ObsDataVector<float> > obserr,
               const eckit::LocalConfiguration &conf);

 private:
  Parameters_ options_;
  eckit::LocalConfiguration unitTestConfig_;

  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override {
    return QCflags::history;
  }

  /// \brief Retrieve all station ids from the ObsAccessor. If string-labelled, ids will be
  /// converted to integers.
  ///
  /// \p stringMap A previously-created map of all unique string station ids across obs spaces,
  /// used as a "central source-of-truth" for the integer values the strings should be mapped to.
  /// \p stationIdVar The parameter used to specify which variable is used to store station ids.
  /// \p obsdb The ObsSpace which station ids will be needed for.
  /// \p obsacc The ObsAccessor used to access observations from the associated ObsSpace.
  std::vector<int> getStationIds(const std::map<std::string, int> &stringMap,
                                 const boost::optional<Variable> &stationIdVar,
                                 const ioda::ObsSpace &obsdb,
                                 const ObsAccessor &obsacc) const;
};


}  // namespace ufo

#endif  // UFO_FILTERS_HISTORYCHECK_H_
