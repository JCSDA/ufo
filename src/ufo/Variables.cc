/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "Variables.h"

#include "oops/base/Variables.h"
#include "eckit/config/LocalConfiguration.h"
#include "Fortran.h"
#include "util/Logger.h"

namespace ufo {

// -----------------------------------------------------------------------------

Variables::Variables(const oops::Variables & oopsvars) {
  oops::Log::debug() << "Variables oops:" << oopsvars << std::endl;
  oops::Log::debug() << "Variables list:" << oopsvars.variables() << std::endl;
  eckit::LocalConfiguration config;
  config.set("nvars", oopsvars.variables().size());
  config.set("variables", oopsvars.variables());
  const eckit::Configuration * conf = &config;
  ufo_var_create_f90(keyVar_, &conf);
  print(oops::Log::debug());
}

// -----------------------------------------------------------------------------

Variables::Variables(const eckit::Configuration & config) {
  oops::Log::debug() << "Variables config:" << config << std::endl;
  const eckit::Configuration * conf = &config;
  ufo_var_create_f90(keyVar_, &conf);
  print(oops::Log::debug());
}

// -----------------------------------------------------------------------------

Variables::Variables(const int keyVar): keyVar_(keyVar) {
}

// -----------------------------------------------------------------------------

Variables::~Variables() {
  ufo_var_delete_f90(keyVar_);
}

// -----------------------------------------------------------------------------

Variables::Variables(const Variables & other) {
  ufo_var_clone_f90(other.keyVar_, keyVar_);
}

// -----------------------------------------------------------------------------

void Variables::print(std::ostream & os) const {
  int nv;
  ufo_var_info_f90(keyVar_, nv);
  os << "ufo::Variables: nvar=" << nv;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
