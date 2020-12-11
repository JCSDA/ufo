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

/// FilterBase: Base class for UFO QC filters

// Filters only need to implement the constructor and the applyFilter method,
// the base class takes care of applying the filter at the pre, prior or post stage.

class FilterBase : public util::Printable {
 public:
  FilterBase(ioda::ObsSpace &, const FilterParametersBase &parameters,
             std::shared_ptr<ioda::ObsDataVector<int> >,
             std::shared_ptr<ioda::ObsDataVector<float> >);
  FilterBase(ioda::ObsSpace &, const eckit::Configuration &,
             std::shared_ptr<ioda::ObsDataVector<int> >,
             std::shared_ptr<ioda::ObsDataVector<float> >);
  ~FilterBase();

  void preProcess();
  void priorFilter(const GeoVaLs &);
  void postFilter(const ioda::ObsVector &, const ObsDiagnostics &);

  oops::Variables requiredVars() const {
    return allvars_.allFromGroup("GeoVaLs").toOopsVariables();}
  oops::Variables requiredHdiagnostics() const {
    return allvars_.allFromGroup("ObsDiag").toOopsVariables();}

 protected:
  ioda::ObsSpace & obsdb_;
  /// For backward compatibility, the full set of filter options (including those required only by
  /// the concrete subclass, not by FilterBase) is stored in this LocalConfiguration object.
  /// It will be removed once all filters have been converted to use Parameters.
  const eckit::LocalConfiguration config_;
  std::shared_ptr<ioda::ObsDataVector<int>> flags_;
  std::shared_ptr<ioda::ObsDataVector<float>> obserr_;
  ufo::Variables allvars_;
  ufo::Variables filtervars_;
  ObsFilterData data_;

 private:
  void doFilter() const;
  virtual void print(std::ostream &) const = 0;
  virtual void applyFilter(const std::vector<bool> &, const Variables &,
                           std::vector<std::vector<bool>> &) const = 0;
  virtual int qcFlag() const = 0;
  bool prior_;
  bool post_;

  // Variables extracted from the filter parameters.
  bool deferToPost_;
  // These Configuration objects will be replaced by appropriate Parameters subclasses later.
  eckit::LocalConfiguration whereConfig_;
  eckit::LocalConfiguration actionConfig_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_FILTERBASE_H_
