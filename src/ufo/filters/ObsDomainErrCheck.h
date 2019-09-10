/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_OBSDOMAINERRCHECK_H_
#define UFO_FILTERS_OBSDOMAINERRCHECK_H_

#include <ostream>
#include <string>
#include <vector>

#include "boost/shared_ptr.hpp"

#include "eckit/config/LocalConfiguration.h"
#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/filters/Variables.h"

namespace eckit {class Configuration;}

namespace ioda {
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

/// Domain check: AMSU-A scattering check and obserr inflation
//  that obs are within domain

// Domain is defined by metadata criteria regardless of obs value.
// If obs value is required, use ObsBoundsCheck.

// The same effect can be achieved with opposite criteria through BlackList,
// the choice is a matter of convenience or which seems more natural.

class ObsDomainErrCheck : public util::Printable,
                       private util::ObjectCounter<ObsDomainErrCheck> {
 public:
  static const std::string classname() {return "ufo::ObsDomainErrCheck";}

  ObsDomainErrCheck(ioda::ObsSpace &, const eckit::Configuration &,
                 boost::shared_ptr<ioda::ObsDataVector<int> >,
                 boost::shared_ptr<ioda::ObsDataVector<float> >);
  ~ObsDomainErrCheck();

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
  const ufo::Variables allvars_;
  const oops::Variables geovars_;
  const oops::Variables diagvars_;
  ioda::ObsDataVector<int> & flags_;
  ioda::ObsDataVector<float> & obserr_;
  float parameter_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_OBSDOMAINERRCHECK_H_
