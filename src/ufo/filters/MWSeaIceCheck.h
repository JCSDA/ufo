/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_MWSEAICECHECK_H_
#define UFO_FILTERS_MWSEAICECHECK_H_

#include <ostream>
#include <string>

#include "boost/shared_ptr.hpp"

#include "eckit/config/LocalConfiguration.h"
#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "ufo/filters/ObsFilterData.h"

namespace ioda {
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

/// MWSeaIceCheck: generic quality control based on observation data only

// Check that observations are within some bounds over some domain

class MWSeaIceCheck : public util::Printable,
                      private util::ObjectCounter<MWSeaIceCheck> {
 public:
  static const std::string classname() {return "ufo::MWSeaIceCheck";}

  MWSeaIceCheck(ioda::ObsSpace &, const eckit::Configuration &,
                boost::shared_ptr<ioda::ObsDataVector<int> >,
                boost::shared_ptr<ioda::ObsDataVector<float> >);
  ~MWSeaIceCheck();

  void preProcess() const {}
  void priorFilter(const GeoVaLs &);
  void postFilter(const ioda::ObsVector &, const ObsDiagnostics &) const {}

  const oops::Variables & requiredGeoVaLs() const {return geovars_;}
  const oops::Variables & requiredHdiagnostics() const {return diagvars_;}

 private:
  void print(std::ostream &) const;

  ioda::ObsSpace & obsdb_;
  ObsFilterData data_;
  const eckit::LocalConfiguration config_;
  const oops::Variables geovars_;
  const oops::Variables diagvars_;
  ioda::ObsDataVector<int> & flags_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_MWSEAICECHECK_H_
