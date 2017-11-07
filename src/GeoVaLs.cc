/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "GeoVaLs.h"

#include "ObsSpace.h"
#include "Fortran.h"
#include "Variables.h"

namespace ufo {
// -----------------------------------------------------------------------------
GeoVaLs::GeoVaLs(const ObsSpace & obsdb, const Variables & var,
                 const util::DateTime & t1, const util::DateTime & t2) {
  const util::DateTime * p1 = &t1;
  const util::DateTime * p2 = &t2;
  qg_obsdb_getgom_f90(obsdb.toFortran(), var.toFortran(), &p1, &p2, keyGVL_);
}
// -----------------------------------------------------------------------------
GeoVaLs::~GeoVaLs() {
  qg_gom_delete_f90(keyGVL_);
}
// -----------------------------------------------------------------------------
void GeoVaLs::zero() {
  qg_gom_zero_f90(keyGVL_);
}
// -----------------------------------------------------------------------------
void GeoVaLs::random() {
  qg_gom_random_f90(keyGVL_);
}
// -----------------------------------------------------------------------------
double GeoVaLs::dot_product_with(const GeoVaLs & other) const {
  double zz;
  qg_gom_dotprod_f90(keyGVL_, other.keyGVL_, zz);
  return zz;
}
// -----------------------------------------------------------------------------
void GeoVaLs::print(std::ostream & os) const {
  os << "GeoVaLs::print not implemented";
}
// -----------------------------------------------------------------------------
}  // namespace ufo
