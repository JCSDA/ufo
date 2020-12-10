/*
 * (C) Copyright 2019 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_TEMPORALTHINNING_H_
#define UFO_FILTERS_TEMPORALTHINNING_H_

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
}

namespace ufo {

class ObsAccessor;
class RecursiveSplitter;
class TemporalThinningParameters;

/// \brief Thin observations so that the retained ones are sufficiently separated in time.
///
/// See TemporalThinningParameters for the documentation of the available parameters.
class TemporalThinning : public FilterBase,
                         private util::ObjectCounter<TemporalThinning> {
 public:
  static const std::string classname() {return "ufo::TemporalThinning";}

  TemporalThinning(ioda::ObsSpace &obsdb, const eckit::Configuration &config,
                   std::shared_ptr<ioda::ObsDataVector<int> > flags,
                   std::shared_ptr<ioda::ObsDataVector<float> > obserr);

  ~TemporalThinning() override;

 private:
  void print(std::ostream &) const override;
  void applyFilter(const std::vector<bool> &, const Variables &,
                   std::vector<std::vector<bool>> &) const override;
  int qcFlag() const override { return QCflags::thinned; }

  ObsAccessor createObsAccessor() const;

  std::vector<bool> identifyThinnedObservations(const std::vector<bool> &apply,
                                                const ObsAccessor &obsAccessor) const;

  boost::optional<std::vector<int>> getObservationPriorities(
      const ObsAccessor &obsAccessor) const;

 private:
  std::unique_ptr<TemporalThinningParameters> options_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_TEMPORALTHINNING_H_
