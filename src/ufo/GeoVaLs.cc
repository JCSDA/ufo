/*
 * (C) Copyright 2017-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/GeoVaLs.h"

#include <iomanip>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/Fortran.h"
#include "ufo/Locations.h"

namespace ufo {

// -----------------------------------------------------------------------------
/*! \brief Constructor given Locations and Variables
 *
 * \details This ufo::GeoVaLs constructor is typically used to initialize
 * GeoVaLs for the full time window (ufo::Locations hold all locations within
 * data assimilation window) and all variables (oops::Variables hold all 
 * variables specified by the ObsOperator as input varialbes. Note that
 * nothing is allocated in the constructor currently, and getValues is
 * responsible for allocation
 *
 */
GeoVaLs::GeoVaLs(const Locations & locs, const oops::Variables & vars)
  : keyGVL_(-1), vars_(vars)
{
  oops::Log::trace() << "GeoVaLs contructor starting" << std::endl;
  const eckit::Configuration * cvar = &vars_.toFortran();
  ufo_geovals_setup_f90(keyGVL_, locs.nobs(), &cvar);
  oops::Log::trace() << "GeoVaLs contructor key = " << keyGVL_ << std::endl;
}

// -----------------------------------------------------------------------------
/*! \brief Constructor for tests
 *
 * \details This ufo::GeoVaLs constructor is typically used in tests, GeoVaLs
 * are read from the file.
 */
GeoVaLs::GeoVaLs(const eckit::Configuration & config,
                 const ioda::ObsSpace & obspace,
                 const oops::Variables & vars)
  : keyGVL_(-1), vars_(vars)
{
  oops::Log::trace() << "GeoVaLs constructor config starting" << std::endl;
  const eckit::Configuration * conf = &config;
  const eckit::Configuration * cvar = &vars_.toFortran();
  ufo_geovals_setup_f90(keyGVL_, 0, &cvar);
  ufo_geovals_read_file_f90(keyGVL_, &conf, obspace, &cvar);
  oops::Log::trace() << "GeoVaLs contructor config key = " << keyGVL_ << std::endl;
}
// -----------------------------------------------------------------------------
/*! \brief Copy constructor */

GeoVaLs::GeoVaLs(const GeoVaLs & other)
  : keyGVL_(-1), vars_(other.vars_)
{
  oops::Log::trace() << "GeoVaLs copy constructor starting" << std::endl;
  const eckit::Configuration * cvar = &vars_.toFortran();
  ufo_geovals_setup_f90(keyGVL_, 0, &cvar);
  ufo_geovals_copy_f90(other.keyGVL_, keyGVL_);
  oops::Log::trace() << "GeoVaLs copy constructor key = " << keyGVL_ << std::endl;
}
// -----------------------------------------------------------------------------
/*? \brief Destructor */
GeoVaLs::~GeoVaLs() {
  oops::Log::trace() << "GeoVaLs destructor starting" << std::endl;
  ufo_geovals_delete_f90(keyGVL_);
  oops::Log::trace() << "GeoVaLs destructor done" << std::endl;
}
// -----------------------------------------------------------------------------
/*! \brief Analytic initialization for GeoVaLs
 *
 * \details This ufo::GeoVaLs constructor was introduced in May, 2018 for use with
 * the interpolation test.   If "analytic_init" is not specified in the
 * configuration then this does nothing.  If "analytic_init" **is** specified, then 
 * the values are replaced by values computed directly from one of several idealized 
 * analytic states.
 *
 * \date May, 2018: Created (M. Miesch, JCSDA)
 * \date June, 2018: Split off from constructor into independent method 
 *                   (M. Miesch, JCSDA)
 */
void GeoVaLs::analytic_init(const Locations & locs,
                            const eckit::Configuration & config)
{
  oops::Log::trace() << "GeoVaLs::analytic_init starting" << std::endl;
  const eckit::Configuration * conf = &config;
  if (config.has("analytic_init")) {
      ufo_geovals_analytic_init_f90(keyGVL_, locs.toFortran(), &conf);
  }
  oops::Log::trace() << "GeoVaLs::analytic_init done" << std::endl;
}
// -----------------------------------------------------------------------------
/*! \brief Zero out the GeoVaLs */
void GeoVaLs::zero() {
  oops::Log::trace() << "GeoVaLs::zero starting" << std::endl;
  ufo_geovals_zero_f90(keyGVL_);
  oops::Log::trace() << "GeoVaLs::zero done" << std::endl;
}
// -----------------------------------------------------------------------------
/*! \brief Absolute value */
void GeoVaLs::abs() {
  oops::Log::trace() << "GeoVaLs::abs starting" << std::endl;
  ufo_geovals_abs_f90(keyGVL_);
  oops::Log::trace() << "GeoVaLs::abs done" << std::endl;
}
// -----------------------------------------------------------------------------
/*! \brief Calculate norm */
double GeoVaLs::norm() const {
  oops::Log::trace() << "GeoVaLs::norm starting" << std::endl;
  double zz;
  ufo_geovals_rms_f90(keyGVL_, zz);
  oops::Log::trace() << "GeoVaLs::norm done" << std::endl;
  return zz;
}
// -----------------------------------------------------------------------------
/*! \brief Randomize GeoVaLs */
void GeoVaLs::random() {
  oops::Log::trace() << "GeoVaLs::random starting" << std::endl;
  ufo_geovals_random_f90(keyGVL_);
  oops::Log::trace() << "GeoVaLs::random done" << std::endl;
}
// -----------------------------------------------------------------------------
/*! \brief Multiply by a constant scalar */
GeoVaLs & GeoVaLs::operator*=(const double zz) {
  oops::Log::trace() << "GeoVaLs::operator*= starting" << std::endl;
  ufo_geovals_scalmult_f90(keyGVL_, zz);
  oops::Log::trace() << "GeoVaLs::operator*= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
/*! \brief Copy operator */
GeoVaLs & GeoVaLs::operator=(const GeoVaLs & rhs) {
  oops::Log::trace() << "GeoVaLs::operator= starting" << std::endl;
  ufo_geovals_assign_f90(keyGVL_, rhs.keyGVL_);
  oops::Log::trace() << "GeoVaLs::operator= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
/*! \brief Add another GeoVaLs */
GeoVaLs & GeoVaLs::operator+=(const GeoVaLs & other) {
  oops::Log::trace() << "GeoVaLs::operator+= starting" << std::endl;
  ufo_geovals_add_f90(keyGVL_, other.keyGVL_);
  oops::Log::trace() << "GeoVaLs::operator+= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
/*! \brief Subtract another GeoVaLs */
GeoVaLs & GeoVaLs::operator-=(const GeoVaLs & other) {
  oops::Log::trace() << "GeoVaLs::operator-= starting" << std::endl;
  ufo_geovals_diff_f90(keyGVL_, other.keyGVL_);
  oops::Log::trace() << "GeoVaLs::operator-= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
/*! \brief Multiply another GeoVaLs */
GeoVaLs & GeoVaLs::operator*=(const GeoVaLs & other) {
  oops::Log::trace() << "GeoVaLs::operator*= starting" << std::endl;
  ufo_geovals_schurmult_f90(keyGVL_, other.keyGVL_);
  oops::Log::trace() << "GeoVaLs::operator*= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
/*! \brief GeoVaLs normalization
 *
 * \details This operator is used to normalize each element of the input GeoVaLs
 * object (LHS) with the rms values of each variable on the RHS, across all
 * locations
 */
GeoVaLs & GeoVaLs::operator/=(const GeoVaLs & other) {
  oops::Log::trace() << "GeoVaLs::operator/= starting" << std::endl;
  ufo_geovals_normalize_f90(keyGVL_, other.keyGVL_);
  oops::Log::trace() << "GeoVaLs::operator/= done" << std::endl;
  return *this;
}
// -----------------------------------------------------------------------------
/*! \brief Scalar product of two GeoVaLs */
double GeoVaLs::dot_product_with(const GeoVaLs & other) const {
  oops::Log::trace() << "GeoVaLs::dot_product_with starting" << std::endl;
  double zz;
  ufo_geovals_dotprod_f90(keyGVL_, other.keyGVL_, zz);
  oops::Log::trace() << "GeoVaLs::dot_product_with done" << std::endl;
  return zz;
}
// -----------------------------------------------------------------------------
/*! \brief Output GeoVaLs to a stream */
void GeoVaLs::print(std::ostream & os) const {
  int nn;
  double zmin, zmax, zrms;

  os << "GeoVaLs: variables = " << vars_ << std::endl;
  for (size_t jv = 0; jv < vars_.size(); ++jv) {
    int nv = jv;
    ufo_geovals_minmaxavg_f90(keyGVL_, nn, nv, zmin, zmax, zrms);
    os << "GeoVaLs: nobs= " << nn << " " << vars_[jv] << " Min=" << zmin << ", Max=" << zmax
       << ", RMS=" << zrms << std::endl;
  }

  /*! Verbose print statement (debug mode)
   *
   * \detail If the min value across all variables is positive, then this may be
   * an error measurement.  If so, compute the rms over the vertical profile and
   * tell the user where the maximum rms value occurs, in terms of the
   * observation number and the variable number.  This is intended to help
   * with debugging.
   */

  if (zmin >= 0.0) {
    double mxval;
    int ivar, iobs;

    ufo_geovals_maxloc_f90(keyGVL_, mxval, iobs, ivar);

    oops::Log::debug() << "GeoVaLs: Maximum Value (vertical rms) = "
                       << std::setprecision(4)
                       << mxval << " for observation = " << iobs
                       << " and variable = " << ivar << std::endl;
  }
}
// -----------------------------------------------------------------------------
/*! \brief Return number of levels for a specified variable */
size_t GeoVaLs::nlevs(const std::string & var) const {
  oops::Log::trace() << "GeoVaLs::nlevs starting" << std::endl;
  int nlevs;
  ufo_geovals_nlevs_f90(keyGVL_, var.size(), var.c_str(), nlevs);
  oops::Log::trace() << "GeoVaLs::nlevs done" << std::endl;
  return nlevs;
}
// -----------------------------------------------------------------------------
/*! \brief Return all values for a specific 2D variable */
void GeoVaLs::get(std::vector<float> & vals, const std::string & var) const {
  oops::Log::trace() << "GeoVaLs::get 2D starting" << std::endl;
  int nlocs;
  ufo_geovals_nlocs_f90(keyGVL_, nlocs);
  ASSERT(vals.size() == nlocs);
  ufo_geovals_get2d_f90(keyGVL_, var.size(), var.c_str(), nlocs, vals[0]);
  oops::Log::trace() << "GeoVaLs::get 2D done" << std::endl;
}
// -----------------------------------------------------------------------------
/*! \brief Return all values for a specific variable and level */
void GeoVaLs::get(std::vector<float> & vals, const std::string & var, const int lev) const {
  oops::Log::trace() << "GeoVaLs::get starting" << std::endl;
  int nlocs;
  ufo_geovals_nlocs_f90(keyGVL_, nlocs);
  ASSERT(vals.size() == nlocs);
  ufo_geovals_get_f90(keyGVL_, var.size(), var.c_str(), lev, nlocs, vals[0]);
  oops::Log::trace() << "GeoVaLs::get done" << std::endl;
}
// -----------------------------------------------------------------------------
/*! \brief Read GeoVaLs from the file */
void GeoVaLs::read(const eckit::Configuration & config,
                   const ioda::ObsSpace & obspace) {
  oops::Log::trace() << "GeoVaLs::read starting" << std::endl;
  const eckit::Configuration * conf = &config;
  const eckit::Configuration * cvar = &vars_.toFortran();
  ufo_geovals_read_file_f90(keyGVL_, &conf, obspace, &cvar);
  oops::Log::trace() << "GeoVaLs::read done" << std::endl;
}
// -----------------------------------------------------------------------------
/*! \brief Write GeoVaLs to the file */
void GeoVaLs::write(const eckit::Configuration & config) const {
  oops::Log::trace() << "GeoVaLs::write starting" << std::endl;
  const eckit::Configuration * conf = &config;
  ufo_geovals_write_file_f90(keyGVL_, &conf);
  oops::Log::trace() << "GeoVaLs::write done" << std::endl;
}
// -----------------------------------------------------------------------------
}  // namespace ufo
