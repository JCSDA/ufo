/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "tools/new_qc/example/Example.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------

Example::Example(ioda::ObsSpace & obsdb, const Parameters_ & parameters,
                 std::shared_ptr<ioda::ObsDataVector<int> > flags,
                 std::shared_ptr<ioda::ObsDataVector<float> >)
  : obsdb_(obsdb), geovars_(), flags_(*flags) {
  oops::Log::trace() << "Example contructor starting" << std::endl;
  ufo_example_create_f90(key_, parameters.toConfiguration(), geovars_);
  oops::Log::debug() << "Example contructor key = " << key_ << std::endl;
}

// -----------------------------------------------------------------------------

Example::~Example() {
  oops::Log::trace() << "Example destructor key = " << key_ << std::endl;
  ufo_example_delete_f90(key_);
}

// -----------------------------------------------------------------------------

void Example::priorFilter(const GeoVaLs & gv) {
  oops::Log::trace() << "Example priorFilter" << std::endl;
  ufo_example_prior_f90(key_, obsdb_, gv.toFortran());
}

// -----------------------------------------------------------------------------

void Example::postFilter(const GeoVaLs &,
                         const ioda::ObsVector & hofxb,
                         const ioda::ObsVector & bias,
                         const ObsDiagnostics & diags) {
  oops::Log::trace() << "Example postFilter" << std::endl;
  ufo_example_post_f90(key_, obsdb_, hofxb.nvars(), hofxb.nlocs(), hofxb.toFortran(),
                       bias.toFortran(), diags.toFortran());
}

// -----------------------------------------------------------------------------

void Example::print(std::ostream & os) const {
  os << "Example::print not yet implemented " << key_;
}
}  // namespace ufo
