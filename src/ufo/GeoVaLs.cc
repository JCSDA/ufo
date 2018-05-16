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
GeoVaLs::GeoVaLs(const Locations & locs, const oops::Variables & vars)
 : keyGVL_(-1), vars_(vars)
{
  oops::Log::trace() << "GeoVaLs contructor starting" << std::endl;
  const eckit::Configuration * cvar = &vars_.toFortran();
  ufo_geovals_setup_f90(keyGVL_, locs.toFortran(), &cvar);
  oops::Log::trace() << "GeoVaLs contructor key = " << keyGVL_ << std::endl;
}
// -----------------------------------------------------------------------------
GeoVaLs::GeoVaLs(const eckit::Configuration & config, const oops::Variables & vars)
 : keyGVL_(-1), vars_(vars)
{
  oops::Log::trace() << "GeoVaLs contructor config starting" << std::endl;
  ufo_geovals_create_f90(keyGVL_);
  int irandom = 0;
  config.get("random", irandom);
  const eckit::Configuration * conf = &config;
  const eckit::Configuration * cvar = &vars_.toFortran();
  if (irandom == 0) {
    ufo_geovals_read_file_f90(keyGVL_, &conf, &cvar);
  } else {
    ufo_geovals_setup_random_f90(keyGVL_, &conf, &cvar);
  }
  oops::Log::trace() << "GeoVaLs contructor config key = " << keyGVL_ << std::endl;
}
// -----------------------------------------------------------------------------
GeoVaLs::GeoVaLs(const Locations & locs, const oops::Variables & vars,
		 const eckit::Configuration & config)
 : keyGVL_(-1), vars_(vars)
{
  oops::Log::trace() << "GeoVaLs contructor config starting" << std::endl;
  ufo_geovals_create_f90(keyGVL_);
  int irandom = 0;
  config.get("random", irandom);
  const eckit::Configuration * conf = &config;
  const eckit::Configuration * cvar = &vars_.toFortran();
  if (irandom == 0) {
    ufo_geovals_read_file_f90(keyGVL_, &conf, &cvar);
  } else {
    ufo_geovals_setup_random_f90(keyGVL_, &conf, &cvar);
  }
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
/*! Absolute value */
void GeoVaLs::abs() {
  ufo_geovals_abs_f90(keyGVL_);
}
// -----------------------------------------------------------------------------
void GeoVaLs::random() {
  ufo_geovals_random_f90(keyGVL_);
}
// -----------------------------------------------------------------------------
GeoVaLs & GeoVaLs::operator*=(const double zz) {
  ufo_geovals_scalmult_f90(keyGVL_, zz);
  return *this;
}
// -----------------------------------------------------------------------------
GeoVaLs & GeoVaLs::operator=(const GeoVaLs & rhs) {
  ufo_geovals_assign_f90(keyGVL_, rhs.keyGVL_);
  return *this;
}
// -----------------------------------------------------------------------------
GeoVaLs & GeoVaLs::operator+=(const GeoVaLs & other) {
  ufo_geovals_add_f90(keyGVL_, other.keyGVL_);
  return *this;
}
// -----------------------------------------------------------------------------
GeoVaLs & GeoVaLs::operator-=(const GeoVaLs & other) {
  ufo_geovals_diff_f90(keyGVL_, other.keyGVL_);
  return *this;
}
// -----------------------------------------------------------------------------
/*! GeoVaLs normalization
 *
 * \details This operator is used to normalize each element of the input GeoVaLs
 * object (LHS) with the rms values of each variable on the RHS, across all 
 * locations
 */
GeoVaLs & GeoVaLs::operator/=(const GeoVaLs & other) {
  ufo_geovals_normalize_f90(keyGVL_, other.keyGVL_);
  return *this;
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
  const eckit::Configuration * cvar = &vars_.toFortran();
  ufo_geovals_read_file_f90(keyGVL_, &conf, &cvar);
}
// -----------------------------------------------------------------------------
void GeoVaLs::write(const eckit::Configuration & config) const {
  const eckit::Configuration * conf = &config;
  ufo_geovals_write_file_f90(keyGVL_, &conf);
}
// -----------------------------------------------------------------------------
}  // namespace ufo
