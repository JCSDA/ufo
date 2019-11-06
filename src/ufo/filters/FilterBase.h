/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_FILTERBASE_H_
#define UFO_FILTERS_FILTERBASE_H_

#include <ostream>
#include <string>
#include <vector>

#include "boost/shared_ptr.hpp"

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
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
  FilterBase(ioda::ObsSpace &, const eckit::Configuration &,
             boost::shared_ptr<ioda::ObsDataVector<int> >,
             boost::shared_ptr<ioda::ObsDataVector<float> >);
  ~FilterBase();

  void preProcess();
  void priorFilter(const GeoVaLs &);
  void postFilter(const ioda::ObsVector &, const ObsDiagnostics &);

  oops::Variables requiredGeoVaLs() const {
    return allvars_.allFromGroup("GeoVaLs").toOopsVariables();}
  oops::Variables requiredHdiagnostics() const {
    return allvars_.allFromGroup("ObsDiag").toOopsVariables();}

 protected:
  ioda::ObsSpace & obsdb_;
  const eckit::LocalConfiguration config_;
  ioda::ObsDataVector<int> & flags_;
  ioda::ObsDataVector<float> & obserr_;
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
};

}  // namespace ufo

#endif  // UFO_FILTERS_FILTERBASE_H_
