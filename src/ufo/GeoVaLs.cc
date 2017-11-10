/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "GeoVaLs.h"

#include "eckit/config/Configuration.h"
#include "ObsSpace.h"
#include "Fortran.h"
#include "Variables.h"
#include "util/Logger.h"

namespace ufo {
// -----------------------------------------------------------------------------
GeoVaLs::GeoVaLs(const ObsSpace & obsdb, const Variables & var,
                 const util::DateTime & t1, const util::DateTime & t2) {
  oops::Log::trace() << "GeoVaLs contructor starting " << t1 << " " << t2 << std::endl;
  const util::DateTime * p1 = &t1;
  const util::DateTime * p2 = &t2;
  ufo_obsdb_getgeovals_f90(obsdb.toFortran(), var.toFortran(), &p1, &p2, keyGVL_);
  oops::Log::trace() << "GeoVaLs contructor key = " << keyGVL_ << std::endl;
}
// -----------------------------------------------------------------------------
GeoVaLs::GeoVaLs(const eckit::Configuration & config, const Variables & var) {
  oops::Log::trace() << "GeoVaLs contructor config starting" << std::endl;
  ufo_geovals_create_f90(keyGVL_);
  const eckit::Configuration * conf = &config;
  ufo_geovals_read_file_f90(keyGVL_, &conf, var.toFortran());
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
void GeoVaLs::read(const eckit::Configuration & config, const Variables & var) {
  const eckit::Configuration * conf = &config;
  ufo_geovals_read_file_f90(keyGVL_, &conf, var.toFortran());
}
// -----------------------------------------------------------------------------
void GeoVaLs::write(const eckit::Configuration & config) const {
  const eckit::Configuration * conf = &config;
  ufo_geovals_write_file_f90(keyGVL_, &conf);
}
// -----------------------------------------------------------------------------
}  // namespace ufo
