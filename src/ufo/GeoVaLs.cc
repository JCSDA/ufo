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

#include "oops/base/Locations.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "oops/util/Range.h"

#include "ufo/GeoVaLs.interface.h"
#include "ufo/ObsTraits.h"
#include "ufo/SampledLocations.h"

namespace ufo {

namespace {

/// \brief Return the Distribution object describing the MPI distribution of locations.
///
/// Technical note: the Distribution is extracted from the SampledLocations objects stored in
/// Locations, so this function throws an exception if `locations` contains no instances of
/// SampledLocations or if not all instances refer to the same Distribution object.
std::shared_ptr<const ioda::Distribution> getSharedDistribution(
    const oops::Locations<ObsTraits> &locations)
{
  if (locations.numSamplingMethods() == 0)
    throw eckit::BadValue("At least one location sampling method must be defined", Here());

  const std::shared_ptr<const ioda::Distribution> dist =
      locations.samplingMethod(0).sampledLocations().distribution();
  for (size_t i = 1; i != locations.numSamplingMethods(); ++i)
    if (locations.samplingMethod(i).sampledLocations().distribution() != dist)
      throw eckit::BadValue("All sampling methods must share the same distribution", Here());
  return dist;
}

}  // namespace

// -----------------------------------------------------------------------------
/*! \brief Deprecated constructor - does not allocate fields.
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

/*! \brief Deprecated constructor given Locations_ and Variables
 *
 * \details Please do not use in any new code. This constructor is currently
 * only used for ObsDiagnostics and will be removed soon. Use the
 *
 *     GeoVaLs(const Locations_ &, const oops::Variables &,
 *             const std::vector<size_t> &)
 *
 * constructor instead.
 * This ufo::GeoVaLs constructor is used to initialize GeoVaLs for specified
 * Locations_ and oops::Variables hold all. Note that nothing is allocated
 * when this constructor is called.
 */
GeoVaLs::GeoVaLs(const Locations_ & locations,
                 const oops::Variables & vars)
  : keyGVL_(-1), vars_(vars), dist_(getSharedDistribution(locations))
{
  oops::Log::trace() << "GeoVaLs contructor starting" << std::endl;

  size_t nlocs;
  std::vector<size_t> numPathsByMethod;
  std::vector<size_t> samplingMethodByVar;
  fillSetupInputs(locations, nlocs, numPathsByMethod, samplingMethodByVar);
  ufo_geovals_partial_setup_f90(keyGVL_, nlocs, vars_, vars_.size(), numPathsByMethod.size(),
                                samplingMethodByVar.data());
  setupSamplingMethods(locations);

  oops::Log::trace() << "GeoVaLs contructor key = " << keyGVL_ << std::endl;
}

// -----------------------------------------------------------------------------
/*! \brief Constructor allocating space for the storage of specified variables.
 *
 * \param locations
 *   Maps variables to sets of paths sampling the observation locations; each variable will be
 *   interpolated along the corresponding set of paths.
 * \param vars
 *   Names of the variables whose values will be stored in the new GeoVaLs object.
 * \param nlevs
 *   Vector whose ith element indicates how many values per interpolation path will be stored
 *   in the GeoVaL corresponding to the ith variable.
 *
 * \details This ufo::GeoVaLs constructor is used in all oops H(x) and DA
 * applications.
 */
GeoVaLs::GeoVaLs(const Locations_ & locations,
                 const oops::Variables & vars, const std::vector<size_t> & nlevs)
  : keyGVL_(-1), vars_(vars), dist_(getSharedDistribution(locations))
{
  oops::Log::trace() << "GeoVaLs constructor starting" << std::endl;

  size_t nlocs;
  std::vector<size_t> numPathsByMethod;
  std::vector<size_t> samplingMethodByVar;
  fillSetupInputs(locations, nlocs, numPathsByMethod, samplingMethodByVar);
  ufo_geovals_setup_f90(keyGVL_, nlocs, vars_, nlevs.size(), nlevs.data(), numPathsByMethod.size(),
                        numPathsByMethod.data(), samplingMethodByVar.data());
  setupSamplingMethods(locations);

  oops::Log::trace() << "GeoVaLs constructor key = " << keyGVL_ << std::endl;
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
  const std::vector<size_t> samplingMethodByVar(vars_.size(), 0);
  ufo_geovals_partial_setup_f90(keyGVL_, 0 /*nlocs*/, vars_, vars_.size(),
                                1 /*num sampling methods*/, samplingMethodByVar.data());
  // only read if there are variables specified
  if (vars_.size() > 0) {
    if (params.filename.value() == boost::none) {
      throw eckit::UserError("geovals requires 'filename' section", Here());
    }
    ufo_geovals_read_file_f90(keyGVL_, params.toConfiguration(), obspace, vars_);
  }
  oops::Log::trace() << "GeoVaLs constructor config key = " << keyGVL_ << std::endl;
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
void GeoVaLs::fillSetupInputs(const Locations_ & locations,
                              size_t &nlocs, std::vector<size_t> & numPathsByMethod,
                              std::vector<size_t> & samplingMethodByVar) const {
  const size_t numSamplingMethods = locations.numSamplingMethods();
  ASSERT(numSamplingMethods > 0);
  nlocs = locations.samplingMethod(0).sampledLocations().nlocs();
  for (size_t i = 1; i < numSamplingMethods; ++i)
    ASSERT(locations.samplingMethod(i).sampledLocations().nlocs() == nlocs);

  numPathsByMethod.resize(numSamplingMethods);
  for (size_t i = 0; i < numSamplingMethods; ++i)
    numPathsByMethod[i] = locations.samplingMethod(i).sampledLocations().size();

  samplingMethodByVar.resize(vars_.size());
  for (size_t i = 0; i < vars_.size(); ++i)
    samplingMethodByVar[i] = locations.samplingMethodIndex(vars_[i]);
}
// -----------------------------------------------------------------------------
void GeoVaLs::setupSamplingMethods(const Locations_ & locations) {
  const size_t nlocs = this->nlocs();
  for (size_t i = 0; i != locations.numSamplingMethods(); ++i) {
    const size_t npaths = locations.samplingMethod(i).sampledLocations().size();
    const std::vector<util::Range<size_t>> &map =
        locations.samplingMethod(i).sampledLocations().pathsGroupedByLocation();
    if (map.empty()) {
      // The location-to-paths map may be omitted only if it's trivial (1-to-1)
      ASSERT(npaths == nlocs);
      ufo_geovals_setup_trivial_sampling_method_f90(keyGVL_, i);
    } else {
      ASSERT(map.size() == nlocs);
      ufo_geovals_setup_sampling_method_f90(keyGVL_, i, npaths, map.size(), map.data());
    }
  }
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
/*! \brief Multiply by a constant scalar for each location */
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
  assert(vars_ == other.vars_);
  auto accumulator = dist_->createAccumulator<double>();
  std::vector<double> this_values, other_values;
  const double missing = util::missingValue(missing);
  const size_t nlocs = this->nlocs();
  std::vector<util::Range<size_t>> profileRangesByLocation;

  // loop over all variables in geovals
  for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
    const size_t nlevs = this->nlevs(vars_[jvar]);
    const size_t nprofiles = this->nprofiles(vars_[jvar]);
    assert(nlevs == other.nlevs(vars_[jvar]));
    assert(nprofiles == other.nprofiles(vars_[jvar]));
    this_values.resize(nprofiles);
    other_values.resize(nprofiles);
    this->getProfileIndicesGroupedByLocation(vars_[jvar], profileRangesByLocation);

    // loop over all levels for this variable
    for (size_t jlev = 0; jlev < nlevs; ++jlev) {
      this->getAtLevel(this_values, vars_[jvar], jlev);
      other.getAtLevel(other_values, vars_[jvar], jlev);
      // loop over all locations
      for (size_t jloc = 0; jloc < nlocs; ++jloc) {
        const util::Range<size_t> &profileRange = profileRangesByLocation[jloc];
        // loop over all paths sampling this location
        double currentLocContribution = 0;
        for (size_t jprofile = profileRange.begin; jprofile < profileRange.end; ++jprofile) {
          if ((this_values[jprofile] != missing) && (other_values[jprofile] != missing)) {
            currentLocContribution += this_values[jprofile] * other_values[jprofile];
          }
        }
        accumulator->addTerm(jloc, currentLocContribution);
      }
    }
  }
  const double dotprod = accumulator->computeResult();
  oops::Log::trace() << "GeoVaLs::dot_product_with done" << std::endl;
  return dotprod;
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
  int nlevs;
  ufo_geovals_nlevs_f90(keyGVL_, var.size(), var.c_str(), nlevs);
  return nlevs;
}
// -----------------------------------------------------------------------------
/*! \brief Return the number of paths along which the specified variable has/will been
 *  interpolated. */
size_t GeoVaLs::nprofiles(const std::string & var) const {
  size_t nprofiles;
  ufo_geovals_nprofiles_f90(keyGVL_, var.size(), var.c_str(), nprofiles);
  return nprofiles;
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
  size_t nprofiles;
  ufo_geovals_nprofiles_f90(keyGVL_, var.size(), var.c_str(), nprofiles);
  ASSERT(vals.size() == nprofiles);
  ufo_geovals_getdouble_f90(keyGVL_, var.size(), var.c_str(), lev, nprofiles, vals[0]);
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
  const size_t nprofiles = this->nprofiles(var);
  ASSERT(vals.size() == nprofiles);
  ufo_geovals_get2d_f90(keyGVL_, var.size(), var.c_str(), nprofiles, vals[0]);
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
void GeoVaLs::getProfile(std::vector<double> & vals,
                         const std::string & var,
                         const int profileIndex) const {
  const size_t nlevs = this->nlevs(var);
  ASSERT(vals.size() == nlevs);
  ASSERT(profileIndex >= 0 && profileIndex < this->nprofiles(var));
  ufo_geovals_get_profile_f90(keyGVL_, var.size(), var.c_str(), profileIndex, nlevs, vals[0]);
}
// -----------------------------------------------------------------------------
void GeoVaLs::getProfile(std::vector<float> & vals,
                         const std::string & var,
                         const int profileIndex) const {
  std::vector <double> doubleVals(vals.size());
  this->getProfile(doubleVals, var, profileIndex);
  this->cast(doubleVals, vals);
}
// -----------------------------------------------------------------------------
void GeoVaLs::getProfile(std::vector<int> & vals,
                         const std::string & var,
                         const int profileIndex) const {
  std::vector <double> doubleVals(vals.size());
  this->getProfile(doubleVals, var, profileIndex);
  this->cast(doubleVals, vals);
}
// -----------------------------------------------------------------------------
/*! \brief Return all values for a specific variable and location */
void GeoVaLs::getAtLocation(std::vector<double> & vals,
                            const std::string & var,
                            const int loc) const {
  ASSERT(this->nprofiles(var) == this->nlocs());
  getProfile(vals, var, loc);
}
// -----------------------------------------------------------------------------
/*! \brief Return all values for a specific variable and location and convert to float */
void GeoVaLs::getAtLocation(std::vector<float> & vals,
                            const std::string & var,
                            const int loc) const {
  ASSERT(this->nprofiles(var) == this->nlocs());
  getProfile(vals, var, loc);
}
// -----------------------------------------------------------------------------
/*! \brief Return all values for a specific variable and location and convert to int */
void GeoVaLs::getAtLocation(std::vector<int> & vals,
                            const std::string & var,
                            const int loc) const {
  ASSERT(this->nprofiles(var) == this->nlocs());
  getProfile(vals, var, loc);
}
// -----------------------------------------------------------------------------
/*! \brief Put double values for a specific variable and level */
void GeoVaLs::putAtLevel(const std::vector<double> & vals,
                         const std::string & var,
                         const int lev) const {
  oops::Log::trace() << "GeoVaLs::putAtLevel(double) starting" << std::endl;
  const size_t np = this->nprofiles(var);
  ASSERT(vals.size() == np);
  ufo_geovals_putdouble_f90(keyGVL_, var.size(), var.c_str(), lev, np, vals[0]);
  oops::Log::trace() << "GeoVaLs::putAtLevel(double) done" << std::endl;
}
// -----------------------------------------------------------------------------
/*! \brief Put float values for a specific variable and level */
void GeoVaLs::putAtLevel(const std::vector<float> & vals,
                         const std::string & var,
                         const int lev) const {
  oops::Log::trace() << "GeoVaLs::putAtLevel(float) starting" << std::endl;
  const size_t nprofiles = this->nprofiles(var);
  ASSERT(vals.size() == nprofiles);
  std::vector<double> doubleVals(vals.begin(), vals.end());
  ufo_geovals_putdouble_f90(keyGVL_, var.size(), var.c_str(), lev, nprofiles, doubleVals[0]);
  oops::Log::trace() << "GeoVaLs::putAtLevel(float) done" << std::endl;
}
// -----------------------------------------------------------------------------
/*! \brief Put int values for a specific variable and level */
void GeoVaLs::putAtLevel(const std::vector<int> & vals,
                         const std::string & var,
                         const int lev) const {
  oops::Log::trace() << "GeoVaLs::putAtLevel(int) starting" << std::endl;
  const size_t nprofiles = this->nprofiles(var);
  ASSERT(vals.size() == nprofiles);
  std::vector<double> doubleVals(vals.begin(), vals.end());
  ufo_geovals_putdouble_f90(keyGVL_, var.size(), var.c_str(), lev, nprofiles, doubleVals[0]);
  oops::Log::trace() << "GeoVaLs::putAtLevel(int) done" << std::endl;
}
// -----------------------------------------------------------------------------
void GeoVaLs::putProfile(const std::vector<double> & vals,
                            const std::string & var,
                            const int profileIndex) const {
  oops::Log::trace() << "GeoVaLs::putProfile(double) starting" << std::endl;
  const size_t nlevs = this->nlevs(var);
  ASSERT(vals.size() == nlevs);
  ASSERT(profileIndex >= 0 && profileIndex < this->nprofiles(var));
  ufo_geovals_put_profile_f90(keyGVL_, var.size(), var.c_str(), profileIndex, nlevs, vals[0]);
  oops::Log::trace() << "GeoVaLs::putProfile(double) done" << std::endl;
}
// -----------------------------------------------------------------------------
void GeoVaLs::putProfile(const std::vector<float> & vals,
                            const std::string & var,
                            const int profileIndex) const {
  oops::Log::trace() << "GeoVaLs::putProfile(float) starting" << std::endl;
  std::vector<double> doubleVals(vals.begin(), vals.end());
  putProfile(doubleVals, var, profileIndex);
  oops::Log::trace() << "GeoVaLs::putProfile(float) done" << std::endl;
}
// -----------------------------------------------------------------------------
void GeoVaLs::putProfile(const std::vector<int> & vals,
                            const std::string & var,
                            const int profileIndex) const {
  oops::Log::trace() << "GeoVaLs::putProfile(int) starting" << std::endl;
  std::vector<double> doubleVals(vals.begin(), vals.end());
  putProfile(doubleVals, var, profileIndex);
  oops::Log::trace() << "GeoVaLs::putProfile(int) done" << std::endl;
}
// -----------------------------------------------------------------------------
/*! \brief Put double values for a specific variable and location */
void GeoVaLs::putAtLocation(const std::vector<double> & vals,
                            const std::string & var,
                            const int loc) const {
  ASSERT(this->nprofiles(var) == this->nlocs());
  putProfile(vals, var, loc);
}
// -----------------------------------------------------------------------------
/*! \brief Put float values for a specific variable and location */
void GeoVaLs::putAtLocation(const std::vector<float> & vals,
                            const std::string & var,
                            const int loc) const {
  ASSERT(this->nprofiles(var) == this->nlocs());
  putProfile(vals, var, loc);
}
// -----------------------------------------------------------------------------
/*! \brief Put int values for a specific variable and location */
void GeoVaLs::putAtLocation(const std::vector<int> & vals,
                            const std::string & var,
                            const int loc) const {
  ASSERT(this->nprofiles(var) == this->nlocs());
  putProfile(vals, var, loc);
}
// -----------------------------------------------------------------------------
void GeoVaLs::getProfileIndicesGroupedByLocation(
    const std::string &var,
    std::vector<util::Range<size_t>> &profileIndicesByLocation) const {
  profileIndicesByLocation.resize(this->nlocs());
  ufo_geovals_get_profile_indices_grouped_by_loc_f90(
        keyGVL_, var.size(), var.c_str(),
        profileIndicesByLocation.size(), profileIndicesByLocation.data());
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
  size_t nlocs;
  ufo_geovals_nlocs_f90(keyGVL_, nlocs);
  return nlocs;
}
// -----------------------------------------------------------------------------
}  // namespace ufo
