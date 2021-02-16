/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_FILTERBASE_H_
#define UFO_FILTERS_FILTERBASE_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "ufo/filters/FilterParametersBase.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/ObsProcessorBase.h"
#include "ufo/filters/Variables.h"

namespace eckit {
  class Configuration;
}

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

/// \brief Base class for UFO QC filters
///
/// Filters only need to implement the constructor and the applyFilter, print and qcFlag methods;
/// the base class takes care of applying the filter at the pre, prior or post stage.

class FilterBase : public ObsProcessorBase {
 public:
  FilterBase(ioda::ObsSpace &, const FilterParametersBaseWithAbstractAction &parameters,
             std::shared_ptr<ioda::ObsDataVector<int> >,
             std::shared_ptr<ioda::ObsDataVector<float> >);
  FilterBase(ioda::ObsSpace &, const eckit::Configuration &,
             std::shared_ptr<ioda::ObsDataVector<int> >,
             std::shared_ptr<ioda::ObsDataVector<float> >);
  ~FilterBase();

 protected:
  /// For backward compatibility, the full set of filter options (including those required only by
  /// the concrete subclass, not by FilterBase) is stored in this LocalConfiguration object.
  /// It will be removed once all filters have been converted to use Parameters.
  const eckit::LocalConfiguration config_;
  ufo::Variables filtervars_;

 private:
  void doFilter() const override;
  void print(std::ostream &) const override = 0;
  virtual void applyFilter(const std::vector<bool> &, const Variables &,
                           std::vector<std::vector<bool>> &) const = 0;
  virtual int qcFlag() const = 0;

  // This Configuration object will be replaced by an appropriate Parameters subclass later.
  eckit::LocalConfiguration whereConfig_;
  std::unique_ptr<FilterActionParametersBase> actionParameters_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_FILTERBASE_H_
