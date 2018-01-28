/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "GeoVaLs.h"

#include "oops/base/Variables.h"
#include "eckit/config/Configuration.h"
#include "Locations.h"
#include "Fortran.h"
#include "util/Logger.h"

namespace ufo {
// -----------------------------------------------------------------------------
GeoVaLs::GeoVaLs(const Locations & locs, const oops::Variables & var) {
  oops::Log::trace() << "GeoVaLs contructor starting" << std::endl;
  const eckit::Configuration * conf = &var.asConfig();
  ufo_geovals_setup_f90(keyGVL_, locs.toFortran(), &conf);
  oops::Log::trace() << "GeoVaLs contructor key = " << keyGVL_ << std::endl;
}
// -----------------------------------------------------------------------------
GeoVaLs::GeoVaLs(const eckit::Configuration & config) {
  oops::Log::trace() << "GeoVaLs contructor config starting" << std::endl;
  ufo_geovals_create_f90(keyGVL_);
  int irandom = 0;
  config.get("random", irandom);
  const eckit::Configuration * conf = &config;
  if (irandom) 
    ufo_geovals_setup_random_f90(keyGVL_, &conf);
  else 
    ufo_geovals_read_file_f90(keyGVL_, &conf);
  oops::Log::trace() << "GeoVaLs contructor config key = " << keyGVL_ << std::endl;
}
// -----------------------------------------------------------------------------
GeoVaLs::~GeoVaLs() {
  ufo_geovals_delete_f90(keyGVL_);
}
// -----------------------------------------------------------------------------
void GeoVaLs::zero() {
  ufo_geovals_zero_f90(keyGVL_);
}
// -----------------------------------------------------------------------------
void GeoVaLs::random() {
  ufo_geovals_random_f90(keyGVL_);
}
// -----------------------------------------------------------------------------
double GeoVaLs::dot_product_with(const GeoVaLs & other) const {
  double zz;
  ufo_geovals_dotprod_f90(keyGVL_, other.keyGVL_, zz);
  return zz;
}
// -----------------------------------------------------------------------------
void GeoVaLs::print(std::ostream & os) const {
  int nn;
  double zmin, zmax, zrms;
  ufo_geovals_minmaxavg_f90(keyGVL_, nn, zmin, zmax, zrms);
  os << "GeoVaLs: nobs= " << nn << " Min=" << zmin << ", Max=" << zmax << ", RMS=" << zrms;
}
// -----------------------------------------------------------------------------
void GeoVaLs::read(const eckit::Configuration & config) {
  const eckit::Configuration * conf = &config;
  ufo_geovals_read_file_f90(keyGVL_, &conf);
}
// -----------------------------------------------------------------------------
void GeoVaLs::write(const eckit::Configuration & config) const {
  const eckit::Configuration * conf = &config;
  ufo_geovals_write_file_f90(keyGVL_, &conf);
}
// -----------------------------------------------------------------------------
}  // namespace ufo
