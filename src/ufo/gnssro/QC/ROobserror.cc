/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/gnssro/QC/ROobserror.h"

#include "eckit/config/Configuration.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/base/ObsFilterBase.h"
#include "oops/interface/ObsFilter.h"
#include "oops/util/Logger.h"
#include "ufo/GeoVaLs.h"
#include "ufo/UfoTrait.h"

namespace ufo {

// -----------------------------------------------------------------------------

ROobserror::ROobserror(ioda::ObsSpace & obsdb,
                       const eckit::Configuration & config,
                       boost::shared_ptr<ioda::ObsDataVector<int> > qc,
                       boost::shared_ptr<ioda::ObsDataVector<float> > oberr)
  : FilterBase(obsdb, config, qc, oberr)
{
  oops::Log::trace() << "ROobserror contructor starting" << std::endl;
  ufo_roobserror_create_f90(key_, obsdb, config);
  oops::Log::trace() << "ROobserror contructor key = " << key_ << std::endl;
}

// -----------------------------------------------------------------------------

ROobserror::~ROobserror() {
  oops::Log::trace() << "ROobserror destructor key = " << key_ << std::endl;
  ufo_roobserror_delete_f90(key_);
}

// -----------------------------------------------------------------------------

void ROobserror::applyFilter(const std::vector<bool> & apply,
                             const Variables & filtervars,
                             std::vector<std::vector<bool>> & flagged) const {
  oops::Log::trace() << "ROobserror using priorFilter" << std::endl;
  flags_->save("FortranQC");    // should pass values to fortran properly
  obserr_->save("FortranERR");  // should pass values to fortran properly
  ufo_roobserror_prior_f90(key_);
  flags_->read("FortranQC");    // should get values from fortran properly
  obserr_->read("FortranERR");  // should get values from fortran properly
}

// -----------------------------------------------------------------------------

void ROobserror::print(std::ostream & os) const {
  os << "ROobserror::print not yet implemented " << key_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
