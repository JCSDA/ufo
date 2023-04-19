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

  /// A string- or integer-valued variable. Observations with different values of that variable will
  /// be thinned separately. If not set all observations will be thinned together.
  oops::OptionalParameter<Variable> categoryVariable{"category_variable", this};
};
}  // namespace ufo

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
}

namespace ufo {
class ObsAccessor;
class RecursiveSplitter;

/// \brief Checks that all surface pressure observations for a given location duting the
/// assimilation window have been calculated using the same source of observed pressure
/// as the reference observation.

class MetOfficePressureConsistencyCheck : public FilterBase,
  private util::ObjectCounter<MetOfficePressureConsistencyCheck> {
 public:
  typedef MetOfficePressureConsistencyCheckParameters Parameters_;

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
  std::vector<size_t>::const_iterator findSeed(
      std::vector<size_t> validObsIds,
      std::vector<util::DateTime> times,
      std::vector<size_t>::const_iterator validObsIndicesBegin,
      std::vector<size_t>::const_iterator validObsIndicesEnd,
      util::DateTime seedTime) const;

 private:
  Parameters_ options_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_METOFFICEPRESSURECONSISTENCYCHECK_H_
