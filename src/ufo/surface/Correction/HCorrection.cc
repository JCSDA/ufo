/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/surface/Correction/HCorrection.h"

#include "eckit/config/Configuration.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/base/Variables.h"
#include "oops/interface/ObsFilter.h"
#include "oops/util/Logger.h"

#include "ufo/ObsDiagnostics.h"
#include "ufo/UfoTrait.h"

namespace ufo {

// -----------------------------------------------------------------------------

HCorrection::HCorrection(const ioda::ObsSpace & obsdb, const eckit::Configuration & config,
                         boost::shared_ptr<ioda::ObsDataVector<int> > flags,
                         boost::shared_ptr<ioda::ObsDataVector<float> >)
  : obsdb_(obsdb), geovars_(), flags_(*flags) {
  oops::Log::trace() << "HCorrection contructor starting" << std::endl;
  const eckit::Configuration * conf = &config;
  ufo_hcorrection_create_f90(key_, conf, geovars_);
  oops::Log::debug() << "HCorrection contructor key = " << key_ << std::endl;
}

// -----------------------------------------------------------------------------

HCorrection::~HCorrection() {
  oops::Log::trace() << "HCorrection destructor key = " << key_ << std::endl;
  ufo_hcorrection_delete_f90(key_);
}

// -----------------------------------------------------------------------------

void HCorrection::priorFilter(const GeoVaLs & gv) const {
  oops::Log::trace() << "HCorrection priorFilter" << std::endl;
  flags_.save("FortranQC");
  ufo_hcorrection_prior_f90(key_, obsdb_, gv.toFortran());
  flags_.read("FortranQC");
}

// -----------------------------------------------------------------------------

void HCorrection::postFilter(const ioda::ObsVector & hofxb, const ObsDiagnostics &) const {
  oops::Log::trace() << "HCorrection postFilter" << std::endl;
  ufo_hcorrection_post_f90(key_, obsdb_, hofxb.nvars(), hofxb.nlocs(), hofxb.toFortran());
}

// -----------------------------------------------------------------------------

void HCorrection::print(std::ostream & os) const {
  os << "HCorrection::print not yet implemented " << key_;
}
}  // namespace ufo
