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
#include "ufo/filters/processWhere.h"
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
///
/// Note: a subclass of FilterBase (called `Filter` from now on) can opt to extract its
/// settings either from a Configuration object or from a subclass of
/// FilterParametersBaseWithAbstractAction (recommended).
///
/// In the former case, `Filter` should provide a constructor with the following signature:
///
///    Filter(ioda::ObsSpace &os, const eckit::Configuration &conf,
///           std::shared_ptr<ioda::ObsDataVector<int> > qcflags,
///           std::shared_ptr<ioda::ObsDataVector<float> > obserrors);
///
/// In the latter case, a subclass of FilterParametersBaseWithAbstractAction holding the settings
/// of the filter in question should be defined. (Usually it can in fact be a subclass of
/// FilterParametersBase; it is necessary to inherit directly from
/// FilterParametersBaseWithAbstractAction only if it is not appropriate for `reject` to be the
/// default action of the filter.) The `Filter` class should then typedef `Parameters_` to the name
/// of that subclass and provide a constructor with the following signature:
///
///    Filter(ioda::ObsSpace &os, const Parameters_ &params,
///           std::shared_ptr<ioda::ObsDataVector<int> > qcflags,
///           std::shared_ptr<ioda::ObsDataVector<float> > obserrors);
///
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
  ufo::Variables filtervars_;

 private:
  void doFilter() const override;
  void print(std::ostream &) const override = 0;
  virtual void applyFilter(const std::vector<bool> &, const Variables &,
                           std::vector<std::vector<bool>> &) const = 0;
  virtual int qcFlag() const = 0;

  std::vector<WhereParameters> whereParameters_;
  std::unique_ptr<FilterActionParametersBase> actionParameters_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_FILTERBASE_H_
