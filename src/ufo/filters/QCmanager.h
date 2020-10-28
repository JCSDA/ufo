/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_FILTERS_QCMANAGER_H_
#define UFO_FILTERS_QCMANAGER_H_

#include <memory>
#include <ostream>

#include "eckit/config/LocalConfiguration.h"
#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/base/Variables.h"
#include "oops/util/Printable.h"

namespace ioda {
  template <typename DATATYPE> class ObsDataVector;
  class ObsVector;
}

namespace ufo {
class GeoVaLs;
class ObsDiagnostics;

class QCmanager : public util::Printable {
 public:
  QCmanager(ioda::ObsSpace &, const eckit::Configuration &,
            std::shared_ptr<ioda::ObsDataVector<int> >,
            std::shared_ptr<ioda::ObsDataVector<float> >);
  ~QCmanager();

  void preProcess() const {}
  void priorFilter(const GeoVaLs &) const {}
  void postFilter(const ioda::ObsVector &, const ObsDiagnostics &) const;

  const oops::Variables & requiredVars() const {return nogeovals_;}
  const oops::Variables & requiredHdiagnostics() const {return nodiags_;}

 private:
  void print(std::ostream &) const;

  ioda::ObsSpace & obsdb_;
  const eckit::LocalConfiguration config_;
  const oops::Variables nogeovals_;
  const oops::Variables nodiags_;
  std::shared_ptr<ioda::ObsDataVector<int>> flags_;
  std::shared_ptr<ioda::ObsDataVector<float>> obserr_;
  const oops::Variables & observed_;
};

}  // namespace ufo

#endif  // UFO_FILTERS_QCMANAGER_H_
