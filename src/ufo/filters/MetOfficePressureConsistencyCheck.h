/*
 * (C) Crown Copyright 2023 Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_METOFFICEPRESSURECONSISTENCYCHECK_H_
#define UFO_FILTERS_METOFFICEPRESSURECONSISTENCYCHECK_H_

#include <memory>
#include <string>
#include <vector>

#include <boost/optional.hpp>

#include "ioda/ObsDataVector.h"
#include "oops/util/ObjectCounter.h"
#include "ufo/filters/FilterBase.h"
#include "ufo/filters/QCflags.h"

namespace ufo {

class MetOfficePressureConsistencyCheckParameters : public FilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(MetOfficePressureConsistencyCheckParameters, FilterParametersBase)

 public:
  /// If set, the filter will use the observation taken as close as possible to \c seed_time, as
  /// the reference observation. If not set, the filter will use the first observation as the
  /// reference observation.
  oops::OptionalParameter<util::DateTime> seedTime{"seed_time", this};
};
}  // namespace ufo

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {
class Variables;
class ObsAccessor;
class RecursiveSplitter;

/// \brief Checks that all surface pressure observations for a given location throughout the
/// assimilation window have been calculated using the same source of observed pressure
/// as the chosen reference observation. The reference observation is chosen based on proximity
/// to a specific time, either the beginning of the window, or one specified by the user.

class MetOfficePressureConsistencyCheck : public FilterBase,
  private util::ObjectCounter<MetOfficePressureConsistencyCheck> {
 public:
  typedef MetOfficePressureConsistencyCheckParameters Parameters_;
  typedef std::vector<size_t>::const_iterator ObsIndexIterator;

  static const std::string classname() {return "ufo::MetOfficePressureConsistencyCheck";}

  MetOfficePressureConsistencyCheck(ioda::ObsSpace &obsdb, const Parameters_ &parameters,
                      std::shared_ptr<ioda::ObsDataVector<int> > flags,
                      std::shared_ptr<ioda::ObsDataVector<float> > obserr);

  ~MetOfficePressureConsistencyCheck() override;

 private:
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  ObsAccessor createObsAccessor() const;
  void print(std::ostream &) const override;
  int qcFlag() const override {return QCflags::black;}

  /// Return an iterator to the reference observation.
  typename MetOfficePressureConsistencyCheck::ObsIndexIterator findReferenceObservation(
      const std::vector<size_t> & validObsIds,
      const std::vector<util::DateTime> & times,
      const std::vector<size_t>::const_iterator & validObsIndicesBegin,
      const std::vector<size_t>::const_iterator & validObsIndicesEnd,
      const util::DateTime & seedTime) const;

 private:
  Parameters_ options_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_METOFFICEPRESSURECONSISTENCYCHECK_H_
