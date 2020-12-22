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
/*! \brief Default constructor - does not allocate fields
*/
GeoVaLs::GeoVaLs(const eckit::mpi::Comm & comm)
  : keyGVL_(-1), comm_(comm)
{
  oops::Log::trace() << "GeoVaLs default constructor starting" << std::endl;
  ufo_geovals_default_constr_f90(keyGVL_);
  oops::Log::trace() << "GeoVaLs default constructor end" << std::endl;
}

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
  : keyGVL_(-1), vars_(vars), comm_(locs.getComm())
{
  oops::Log::trace() << "GeoVaLs contructor starting" << std::endl;
  ufo_geovals_setup_f90(keyGVL_, locs.nobs(), vars_);
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
  : keyGVL_(-1), vars_(vars), comm_(obspace.comm())
{
  oops::Log::trace() << "GeoVaLs constructor config starting" << std::endl;
  ufo_geovals_setup_f90(keyGVL_, 0, vars_);
  // only read if there are variables specified
  if (vars_.size() > 0)  ufo_geovals_read_file_f90(keyGVL_, config, obspace, vars_);
  oops::Log::trace() << "GeoVaLs contructor config key = " << keyGVL_ << std::endl;
}
// -----------------------------------------------------------------------------
/*! \brief Construct a new GeoVaLs with just one location
*
* \details This ufo::GeoVaLs constructor takes a GeoVaLs object and an index to
* create a new GeoVaLs with just one location
*/
GeoVaLs::GeoVaLs(const GeoVaLs & other, const int & index)
  : keyGVL_(-1), vars_(other.vars_), comm_(other.comm_)
{
  oops::Log::trace() << "GeoVaLs copy one GeoVaLs constructor starting" << std::endl;
  ufo_geovals_setup_f90(keyGVL_, 1, vars_);
  int fort_index = index + 1;  // Fortran numbers from 1
  ufo_geovals_copy_one_f90(keyGVL_, other.keyGVL_, fort_index);
  oops::Log::trace() << "GeoVaLs copy one GeoVaLs constructor key = " << keyGVL_ << std::endl;
}
// -----------------------------------------------------------------------------
/*! \brief Copy constructor */

GeoVaLs::GeoVaLs(const GeoVaLs & other)
  : keyGVL_(-1), vars_(other.vars_), comm_(other.comm_)
{
  oops::Log::trace() << "GeoVaLs copy constructor starting" << std::endl;
  ufo_geovals_setup_f90(keyGVL_, 0, vars_);
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
/*! \brief Zero out the GeoVaLs */
void GeoVaLs::zero() {
  oops::Log::trace() << "GeoVaLs::zero starting" << std::endl;
  ufo_geovals_zero_f90(keyGVL_);
  oops::Log::trace() << "GeoVaLs::zero done" << std::endl;
}
// -----------------------------------------------------------------------------
/*! \brief Reorder GeoVaLs in vertical dimension based on vertical coordinate variable */
void GeoVaLs::reorderzdir(const std::string & varname, const std::string & vardir) {
  oops::Log::trace() << "GeoVaLs::reorderzdir starting" << std::endl;
  ufo_geovals_reorderzdir_f90(keyGVL_, varname.size(), varname.c_str(),
                              vardir.size(), vardir.c_str());
  oops::Log::trace() << "GeoVaLs::reorderzdir done" << std::endl;
}
// -----------------------------------------------------------------------------
/*! \brief Calculate rms */
double GeoVaLs::rms() const {
  oops::Log::trace() << "GeoVaLs::rms starting" << std::endl;
  double zz;
  ufo_geovals_rms_f90(keyGVL_, zz);
  oops::Log::trace() << "GeoVaLs::rms done" << std::endl;
  return zz;
}
// -----------------------------------------------------------------------------
/*! \brief Calculate normalized rms */
double GeoVaLs::normalizedrms(const GeoVaLs & other) const {
  oops::Log::trace() << "GeoVaLs::normalizerms starting" << std::endl;
  GeoVaLs temp_gval(*this);
  ufo_geovals_normalize_f90(temp_gval.keyGVL_, other.keyGVL_);
  double zz = temp_gval.rms();
  oops::Log::trace() << "GeoVaLs::normalizerms done" << std::endl;
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
/*! \brief Multiply by a constant scalar for each profile */
GeoVaLs & GeoVaLs::operator*=(const std::vector<float> & vals) {
  oops::Log::trace() << "GeoVaLs::operator*= starting" << std::endl;
  size_t nlocs;
  ufo_geovals_nlocs_f90(keyGVL_, nlocs);
  oops::Log::trace() << "vals, nlocs = " << vals.size() << "   " << nlocs << std::endl;
  ASSERT(vals.size() == nlocs);
  ufo_geovals_profmult_f90(keyGVL_, nlocs, vals[0]);
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
/*! \brief Scalar product of two GeoVaLs */
double GeoVaLs::dot_product_with(const GeoVaLs & other) const {
  oops::Log::trace() << "GeoVaLs::dot_product_with starting" << std::endl;
  double zz;
  ufo_geovals_dotprod_f90(keyGVL_, other.keyGVL_, zz, comm_);
  oops::Log::trace() << "GeoVaLs::dot_product_with done" << std::endl;
  return zz;
}
// -----------------------------------------------------------------------------
/*! \brief Split two GeoVaLs */
void GeoVaLs::split(GeoVaLs & other1, GeoVaLs & other2) const {
  oops::Log::trace() << "GeoVaLs::split GeoVaLs into 2" << std::endl;
  ufo_geovals_split_f90(keyGVL_, other1.keyGVL_, other2.keyGVL_);
  oops::Log::trace() << "GeoVaLs::split GeoVaLs into 2" << std::endl;
  return;
}
// -----------------------------------------------------------------------------
/*! \brief Merge two GeoVaLs */
void GeoVaLs::merge(const GeoVaLs & other1, const GeoVaLs & other2) {
  oops::Log::trace() << "GeoVaLs::merge 2 GeoVaLs" << std::endl;
  ufo_geovals_merge_f90(keyGVL_, other1.keyGVL_, other2.keyGVL_);
  oops::Log::trace() << "GeoVaLs::merge 2 GeoVaLs" << std::endl;
  return;
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
  size_t nlocs;
  ufo_geovals_nlocs_f90(keyGVL_, nlocs);
  ASSERT(vals.size() == nlocs);
  ufo_geovals_get2d_f90(keyGVL_, var.size(), var.c_str(), nlocs, vals[0]);
  oops::Log::trace() << "GeoVaLs::get 2D done" << std::endl;
}
// -----------------------------------------------------------------------------
/*! \brief Return all values for a specific variable and level */
void GeoVaLs::get(std::vector<float> & vals, const std::string & var, const int lev) const {
  oops::Log::trace() << "GeoVaLs::get starting" << std::endl;
  size_t nlocs;
  ufo_geovals_nlocs_f90(keyGVL_, nlocs);
  ASSERT(vals.size() == nlocs);
  ufo_geovals_get_f90(keyGVL_, var.size(), var.c_str(), lev, nlocs, vals[0]);
  oops::Log::trace() << "GeoVaLs::get done" << std::endl;
}
// -----------------------------------------------------------------------------
/*! \brief Return all values for a specific variable and level */
void GeoVaLs::get(std::vector<double> & vals, const std::string & var, const int lev) const {
  oops::Log::trace() << "GeoVaLs::get starting" << std::endl;
  size_t nlocs;
  ufo_geovals_nlocs_f90(keyGVL_, nlocs);
  ASSERT(vals.size() == nlocs);
  ufo_geovals_getdouble_f90(keyGVL_, var.size(), var.c_str(), lev, nlocs, vals[0]);
  oops::Log::trace() << "GeoVaLs::get done" << std::endl;
}
// -----------------------------------------------------------------------------
/*! \brief Put values for a specific variable and level */
void GeoVaLs::put(const std::vector<double> & vals, const std::string & var, const int lev) const {
  oops::Log::trace() << "GeoVaLs::put starting" << std::endl;
  size_t nlocs;
  ufo_geovals_nlocs_f90(keyGVL_, nlocs);
  ASSERT(vals.size() == nlocs);
  ufo_geovals_putdouble_f90(keyGVL_, var.size(), var.c_str(), lev, nlocs, vals[0]);
  oops::Log::trace() << "GeoVaLs::get done" << std::endl;
}
// -----------------------------------------------------------------------------
/*! \brief Read GeoVaLs from the file */
void GeoVaLs::read(const eckit::Configuration & config,
                   const ioda::ObsSpace & obspace) {
  oops::Log::trace() << "GeoVaLs::read starting" << std::endl;
  ufo_geovals_read_file_f90(keyGVL_, config, obspace, vars_);
  oops::Log::trace() << "GeoVaLs::read done" << std::endl;
}
// -----------------------------------------------------------------------------
/*! \brief Write GeoVaLs to the file */
void GeoVaLs::write(const eckit::Configuration & config) const {
  oops::Log::trace() << "GeoVaLs::write starting" << std::endl;
  ufo_geovals_write_file_f90(keyGVL_, config, comm_);
  oops::Log::trace() << "GeoVaLs::write done" << std::endl;
}
// -----------------------------------------------------------------------------
/*! \brief Return the number of geovals */
size_t GeoVaLs::nlocs() const {
  oops::Log::trace() << "GeoVaLs::nlocs starting" << std::endl;
  size_t nlocs;
  ufo_geovals_nlocs_f90(keyGVL_, nlocs);
  oops::Log::trace() << "GeoVaLs::nlocs done" << std::endl;
  return nlocs;
}
// -----------------------------------------------------------------------------
}  // namespace ufo
