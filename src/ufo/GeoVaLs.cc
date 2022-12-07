/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/GeoVaLs.h"

#include <cassert>
#include <iomanip>
#include <utility>
#include <vector>

#include "eckit/exception/Exceptions.h"

#include "ioda/distribution/Accumulator.h"
#include "ioda/distribution/Distribution.h"
#include "ioda/ObsSpace.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.interface.h"
#include "ufo/Locations.h"

namespace ufo {

// -----------------------------------------------------------------------------
/*! \brief Deprecated default constructor - does not allocate fields.
 *
 * \details Please do not use this constructor in new code.
 */
GeoVaLs::GeoVaLs(std::shared_ptr<const ioda::Distribution> dist,
                 const oops::Variables & vars)
  : keyGVL_(-1), vars_(vars), dist_(std::move(dist))
{
  oops::Log::trace() << "GeoVaLs default constructor starting" << std::endl;
  ufo_geovals_default_constr_f90(keyGVL_);
  oops::Log::trace() << "GeoVaLs default constructor end" << std::endl;
}

/*! \brief Deprecated constructor given Locations and Variables
 *
 * \details Please do not use in any new code. This constructor is currently
 * only used for ObsDiagnostics and will be removed soon. Use the
 * GeoVaLs(const Locations &, const oops::Variables &, const std::vector<size_t> &)
 * constructor instead.
 * This ufo::GeoVaLs constructor is used to initialize GeoVaLs for specified
 * ufo::Locations and oops::Variables hold all. Note that nothing is allocated
 * when this constructor is called.
 */
GeoVaLs::GeoVaLs(const Locations & locs, const oops::Variables & vars)
  : keyGVL_(-1), vars_(vars), dist_(locs.distribution())
{
  oops::Log::trace() << "GeoVaLs contructor starting" << std::endl;
  ufo_geovals_partial_setup_f90(keyGVL_, locs.size(), vars_);
  oops::Log::trace() << "GeoVaLs contructor key = " << keyGVL_ << std::endl;
}

// -----------------------------------------------------------------------------
/*! \brief Allocating constructor for specified Locations \p locs, Variables \p vars
 * and number of levels \p nlevs
 *
 * \details This ufo::GeoVaLs constructor is used in all oops H(x) and DA
 * applications.
 * Sizes of GeoVaLs for i-th variable at a single location are defined by i-th value
 * of \p nlevs.
 */
GeoVaLs::GeoVaLs(const Locations & locs, const oops::Variables & vars,
                 const std::vector<size_t> & nlevs)
  : keyGVL_(-1), vars_(vars), dist_(locs.distribution())
{
  oops::Log::trace() << "GeoVaLs contructor starting" << std::endl;
  ufo_geovals_setup_f90(keyGVL_, locs.size(), vars_, nlevs.size(), nlevs[0]);
  oops::Log::trace() << "GeoVaLs contructor key = " << keyGVL_ << std::endl;
}

// -----------------------------------------------------------------------------
/*! \brief Constructor for tests
 *
 * \details This ufo::GeoVaLs constructor is typically used in tests, GeoVaLs
 * are read from the file.
 */
GeoVaLs::GeoVaLs(const Parameters_ & params,
                 const ioda::ObsSpace & obspace,
                 const oops::Variables & vars)
  : keyGVL_(-1), vars_(vars), dist_(obspace.distribution())
{
  oops::Log::trace() << "GeoVaLs constructor config starting" << std::endl;
  ufo_geovals_partial_setup_f90(keyGVL_, 0, vars_);
  // only read if there are variables specified
  if (vars_.size() > 0) {
    if (params.filename.value() == boost::none) {
      throw eckit::UserError("geovals requires 'filename' section", Here());
    }
    ufo_geovals_read_file_f90(keyGVL_, params.toConfiguration(), obspace, vars_);
  }
  oops::Log::trace() << "GeoVaLs contructor config key = " << keyGVL_ << std::endl;
}
// -----------------------------------------------------------------------------
/*! \brief Construct a new GeoVaLs with just one location
*
* \details This ufo::GeoVaLs constructor takes a GeoVaLs object and an index to
* create a new GeoVaLs with just one location
*/
GeoVaLs::GeoVaLs(const GeoVaLs & other, const int & index)
  : keyGVL_(-1), vars_(other.vars_), dist_(other.dist_)
{
  oops::Log::trace() << "GeoVaLs copy one GeoVaLs constructor starting" << std::endl;
  ufo_geovals_copy_one_f90(keyGVL_, other.keyGVL_, index);
  oops::Log::trace() << "GeoVaLs copy one GeoVaLs constructor key = " << keyGVL_ << std::endl;
}
// -----------------------------------------------------------------------------
/*! \brief Copy constructor */

GeoVaLs::GeoVaLs(const GeoVaLs & other)
  : keyGVL_(-1), vars_(other.vars_), dist_(other.dist_)
{
  oops::Log::trace() << "GeoVaLs copy constructor starting" << std::endl;
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
void GeoVaLs::allocate(const int & nlevels, const oops::Variables & vars)
{
  oops::Log::trace() << "GeoVaLs::allocate starting" << std::endl;
  ufo_geovals_allocate_f90(keyGVL_, nlevels, vars);
  oops::Log::trace() << "GeoVaLs::allocate done" << std::endl;
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
  const size_t nlocs = this->nlocs();
  assert(nlocs == other.nlocs());
  assert(vars_ == other.vars_);
  auto accumulator = dist_->createAccumulator<double>();
  std::vector<double> this_values(nlocs), other_values(nlocs);
  const double missing = util::missingValue(missing);
  // loop over all variables in geovals
  for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
    const size_t nlevs = this->nlevs(vars_[jvar]);
    assert(nlevs == other.nlevs(vars_[jvar]));
    // loop over all levels for this variable
    for (size_t jlev = 0; jlev < nlevs; ++jlev) {
      this->getAtLevel(this_values, vars_[jvar], jlev);
      other.getAtLevel(other_values, vars_[jvar], jlev);
      // loop over all locations
      for (size_t jloc = 0; jloc < nlocs; ++jloc) {
        if ((this_values[jloc] != missing) && (other_values[jloc] != missing)) {
          accumulator->addTerm(jloc, this_values[jloc]*other_values[jloc]);
        }
      }
    }
  }
  const double dotprod = accumulator->computeResult();
  oops::Log::trace() << "GeoVaLs::dot_product_with done" << std::endl;
  return dotprod;
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
  /// Call method to get double values (Fortran data structure stores data in double)
  /// and convert to floats
  std::vector<double> doubleVals(vals.size());
  this->get(doubleVals, var);
  this->cast(doubleVals, vals);
  oops::Log::trace() << "GeoVaLs::get 2D(float) done" << std::endl;
}
// -----------------------------------------------------------------------------
/*! \brief Return all values for a specific variable and level */
void GeoVaLs::getAtLevel(std::vector<double> & vals, const std::string & var, const int lev) const {
  oops::Log::trace() << "GeoVaLs::getAtLevel(double) starting" << std::endl;
  size_t nlocs;
  ufo_geovals_nlocs_f90(keyGVL_, nlocs);
  ASSERT(vals.size() == nlocs);
  ufo_geovals_getdouble_f90(keyGVL_, var.size(), var.c_str(), lev, nlocs, vals[0]);
  oops::Log::trace() << "GeoVaLs::getAtLevel(double) done" << std::endl;
}
// -----------------------------------------------------------------------------
/*! \brief Return all values for a specific variable and level and convert to float */
void GeoVaLs::getAtLevel(std::vector<float> & vals, const std::string & var, const int lev) const {
  oops::Log::trace() << "GeoVaLs::getAtLevel(float) starting" << std::endl;
  std::vector<double> doubleVals(vals.size());
  this->getAtLevel(doubleVals, var, lev);
  this->cast(doubleVals, vals);
  oops::Log::trace() << "GeoVaLs::getAtLevel(float) done" << std::endl;
}
// -----------------------------------------------------------------------------
/*! \brief Return all values for a specific variable and level and convert to int */
void GeoVaLs::getAtLevel(std::vector<int> & vals, const std::string & var, const int lev) const {
  oops::Log::trace() << "GeoVaLs::getAtLevel(int) starting" << std::endl;
  std::vector<double> doubleVals(vals.size());
  this->getAtLevel(doubleVals, var, lev);
  this->cast(doubleVals, vals);
  oops::Log::trace() << "GeoVaLs::getAtLevel(int) done" << std::endl;
}
// -----------------------------------------------------------------------------
/*! \brief Return all values for a specific 2D variable */
void GeoVaLs::get(std::vector<double> & vals, const std::string & var) const {
  oops::Log::trace() << "GeoVaLs::get 2D starting" << std::endl;
  size_t nlocs;
  ufo_geovals_nlocs_f90(keyGVL_, nlocs);
  ASSERT(vals.size() == nlocs);
  ufo_geovals_get2d_f90(keyGVL_, var.size(), var.c_str(), nlocs, vals[0]);
  oops::Log::trace() << "GeoVaLs::get 2D(double) done" << std::endl;
}
// -----------------------------------------------------------------------------
/*! \brief Return all values for a specific 2D variable */
void GeoVaLs::get(std::vector<int> & vals, const std::string & var) const {
  oops::Log::trace() << "GeoVaLs::get 2D starting" << std::endl;
  /// Call method to get double values (Fortran data structure stores data in double)
  /// and convert to ints
  std::vector<double> doubleVals(vals.size());
  this->get(doubleVals, var);
  this->cast(doubleVals, vals);
  oops::Log::trace() << "GeoVaLs::get 2D(int) done" << std::endl;
}
// -----------------------------------------------------------------------------
/*! \brief Return all values for a specific variable and location */
void GeoVaLs::getAtLocation(std::vector<double> & vals,
                            const std::string & var,
                            const int loc) const {
  oops::Log::trace() << "GeoVaLs::getAtLocation(double) starting" << std::endl;
  const size_t nlevs = this->nlevs(var);
  ASSERT(vals.size() == nlevs);
  ASSERT(loc >= 0 && loc < this->nlocs());
  ufo_geovals_get_loc_f90(keyGVL_, var.size(), var.c_str(), loc, nlevs, vals[0]);
  oops::Log::trace() << "GeoVaLs::getAtLocation(double) done" << std::endl;
}
// -----------------------------------------------------------------------------
/*! \brief Return all values for a specific variable and location and convert to float */
void GeoVaLs::getAtLocation(std::vector<float> & vals,
                            const std::string & var,
                            const int loc) const {
  oops::Log::trace() << "GeoVaLs::getAtLocation(float) starting" << std::endl;
  std::vector <double> doubleVals(vals.size());
  this->getAtLocation(doubleVals, var, loc);
  this->cast(doubleVals, vals);
  oops::Log::trace() << "GeoVaLs::getAtLocation(float) done" << std::endl;
}
// -----------------------------------------------------------------------------
/*! \brief Return all values for a specific variable and location and convert to int */
void GeoVaLs::getAtLocation(std::vector<int> & vals,
                            const std::string & var,
                            const int loc) const {
  oops::Log::trace() << "GeoVaLs::getAtLocation(int) starting" << std::endl;
  std::vector <double> doubleVals(vals.size());
  this->getAtLocation(doubleVals, var, loc);
  this->cast(doubleVals, vals);
  oops::Log::trace() << "GeoVaLs::getAtLocation(int) done" << std::endl;
}
// -----------------------------------------------------------------------------
/*! \brief Put double values for a specific variable and level */
void GeoVaLs::putAtLevel(const std::vector<double> & vals,
                         const std::string & var,
                         const int lev) const {
  oops::Log::trace() << "GeoVaLs::putAtLevel(double) starting" << std::endl;
  size_t nlocs;
  ufo_geovals_nlocs_f90(keyGVL_, nlocs);
  ASSERT(vals.size() == nlocs);
  ufo_geovals_putdouble_f90(keyGVL_, var.size(), var.c_str(), lev, nlocs, vals[0]);
  oops::Log::trace() << "GeoVaLs::putAtLevel(double) done" << std::endl;
}
// -----------------------------------------------------------------------------
/*! \brief Put float values for a specific variable and level */
void GeoVaLs::putAtLevel(const std::vector<float> & vals,
                         const std::string & var,
                         const int lev) const {
  oops::Log::trace() << "GeoVaLs::putAtLevel(float) starting" << std::endl;
  size_t nlocs;
  ufo_geovals_nlocs_f90(keyGVL_, nlocs);
  ASSERT(vals.size() == nlocs);
  std::vector<double> doubleVals(vals.begin(), vals.end());
  ufo_geovals_putdouble_f90(keyGVL_, var.size(), var.c_str(), lev, nlocs, doubleVals[0]);
  oops::Log::trace() << "GeoVaLs::putAtLevel(float) done" << std::endl;
}
// -----------------------------------------------------------------------------
/*! \brief Put int values for a specific variable and level */
void GeoVaLs::putAtLevel(const std::vector<int> & vals,
                         const std::string & var,
                         const int lev) const {
  oops::Log::trace() << "GeoVaLs::putAtLevel(int) starting" << std::endl;
  size_t nlocs;
  ufo_geovals_nlocs_f90(keyGVL_, nlocs);
  ASSERT(vals.size() == nlocs);
  std::vector<double> doubleVals(vals.begin(), vals.end());
  ufo_geovals_putdouble_f90(keyGVL_, var.size(), var.c_str(), lev, nlocs, doubleVals[0]);
  oops::Log::trace() << "GeoVaLs::putAtLevel(int) done" << std::endl;
}
/*! \brief Put double values for a specific variable and location */
void GeoVaLs::putAtLocation(const std::vector<double> & vals,
                            const std::string & var,
                            const int loc) const {
  oops::Log::trace() << "GeoVaLs::putAtLocation(double) starting" << std::endl;
  const size_t nlevs = this->nlevs(var);
  ASSERT(vals.size() == nlevs);
  ASSERT(loc >= 0 && loc < this->nlocs());
  ufo_geovals_put_loc_f90(keyGVL_, var.size(), var.c_str(), loc, nlevs, vals[0]);
  oops::Log::trace() << "GeoVaLs::putAtLocation(double) done" << std::endl;
}
/*! \brief Put float values for a specific variable and location */
void GeoVaLs::putAtLocation(const std::vector<float> & vals,
                            const std::string & var,
                            const int loc) const {
  oops::Log::trace() << "GeoVaLs::putAtLocation(float) starting" << std::endl;
  const size_t nlevs = this->nlevs(var);
  ASSERT(vals.size() == nlevs);
  ASSERT(loc >= 0 && loc < this->nlocs());
  std::vector<double> doubleVals(vals.begin(), vals.end());
  ufo_geovals_put_loc_f90(keyGVL_, var.size(), var.c_str(), loc, nlevs, doubleVals[0]);
  oops::Log::trace() << "GeoVaLs::putAtLocation(float) done" << std::endl;
}
/*! \brief Put int values for a specific variable and location */
void GeoVaLs::putAtLocation(const std::vector<int> & vals,
                            const std::string & var,
                            const int loc) const {
  oops::Log::trace() << "GeoVaLs::putAtLocation(int) starting" << std::endl;
  const size_t nlevs = this->nlevs(var);
  ASSERT(vals.size() == nlevs);
  ASSERT(loc >= 0 && loc < this->nlocs());
  std::vector<double> doubleVals(vals.begin(), vals.end());
  ufo_geovals_put_loc_f90(keyGVL_, var.size(), var.c_str(), loc, nlevs, doubleVals[0]);
  oops::Log::trace() << "GeoVaLs::putAtLocation(int) done" << std::endl;
}
// -----------------------------------------------------------------------------
void GeoVaLs::fill(const std::string &name, const ConstVectorRef<size_t> &indx,
                   const ConstMatrixRef<double> &vals, const bool levelsTopDown) {
  oops::Log::trace() << "GeoVaLs::fill starting" << std::endl;
  const size_t npts = indx.size();
  const size_t nlev = vals.cols();
  std::vector<int> findx(indx.size());
  for (Eigen::Index jj = 0; jj < indx.size(); ++jj) findx[jj] = indx[jj];

  ufo_geovals_fill_f90(keyGVL_, name.size(), name.data(),
                       npts, findx.data(), nlev, vals.data(), levelsTopDown);

  oops::Log::trace() << "GeoVaLs::fill done" << std::endl;
}
// -----------------------------------------------------------------------------
void GeoVaLs::fillAD(const std::string &name, const ConstVectorRef<size_t> &indx,
                     MatrixRef<double> vals, const bool levelsTopDown) const {
  oops::Log::trace() << "GeoVaLs::fillAD starting" << std::endl;
  const size_t npts = indx.size();
  const size_t nlev = vals.cols();
  std::vector<int> findx(indx.size());
  for (Eigen::Index jj = 0; jj < indx.size(); ++jj) findx[jj] = indx[jj];

  ufo_geovals_fillad_f90(keyGVL_, name.size(), name.data(),
                            npts, findx.data(), nlev, vals.data(), levelsTopDown);

  oops::Log::trace() << "GeoVaLs::fillAD done" << std::endl;
}
// -----------------------------------------------------------------------------
/*! \brief Read GeoVaLs from the file */
void GeoVaLs::read(const Parameters_ & params,
                   const ioda::ObsSpace & obspace) {
  oops::Log::trace() << "GeoVaLs::read starting" << std::endl;
  if (params.filename.value() == boost::none) {
    throw eckit::UserError("geovals requires 'filename' section", Here());
  }
  ufo_geovals_read_file_f90(keyGVL_, params.toConfiguration(), obspace, vars_);
  oops::Log::trace() << "GeoVaLs::read done" << std::endl;
}
// -----------------------------------------------------------------------------
/*! \brief Write GeoVaLs to the file */
void GeoVaLs::write(const Parameters_ & params) const {
  oops::Log::trace() << "GeoVaLs::write starting" << std::endl;
  if (params.filename.value() == boost::none) {
    throw eckit::UserError("geovals requires 'filename' section", Here());
  }
  ufo_geovals_write_file_f90(keyGVL_, params.toConfiguration(), dist_->rank());
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
