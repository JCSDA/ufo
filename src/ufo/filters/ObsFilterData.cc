/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/ObsFilterData.h"

#include <string>
#include <vector>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/obsfunctions/ObsFunction.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------
ObsFilterData::ObsFilterData(ioda::ObsSpace & obsdb)
  : obsdb_(obsdb), gvals_(NULL), ovecs_(), diags_(NULL) {
  oops::Log::trace() << "ObsFilterData created" << std::endl;
}

// -----------------------------------------------------------------------------
ObsFilterData::~ObsFilterData() {
  oops::Log::trace() << "ObsFilterData destructed" << std::endl;
}

// -----------------------------------------------------------------------------
/*! Associates GeoVaLs with this ObsFilterData (after this call GeoVaLs are available) */
void ObsFilterData::associate(const GeoVaLs & gvals) {
  gvals_ = &gvals;
}

// -----------------------------------------------------------------------------
/*! Associates H(x) ObsVector with this ObsFilterData */
void ObsFilterData::associate(const ioda::ObsVector & hofx, const std::string & name) {
  ovecs_[name] = &hofx;
}

// -----------------------------------------------------------------------------
/*! Associates ObsDiagnostics coming from ObsOperator with this ObsFilterData */
void ObsFilterData::associate(const ObsDiagnostics & diags) {
  diags_ = &diags;
}

// -----------------------------------------------------------------------------
/*! Returns number of observation locations */
size_t ObsFilterData::nlocs() const {
  return obsdb_.nlocs();
}

// -----------------------------------------------------------------------------
/*! Checks if requested data exists in ObsFilterData
 *  \param varname is a name of a variable requested
 *  \return true if the variable is available for access from ObsFilterData,
 *          false otherwise
 */
bool ObsFilterData::has(const Variable & varname) const {
  const std::string var = varname.variable();
  const std::string grp = varname.group();
  if (grp == "GeoVaLs") {
    return (gvals_ && gvals_->has(var));
  } else if (grp == "ObsFunction") {
    return ObsFunctionFactory::functionExists(var);
  } else if (grp == "ObsDiag") {
    return (diags_ && diags_->has(var));
  } else {
    return this->hasVector(grp, var) || obsdb_.has(grp, var);
  }
  return false;
}

// -----------------------------------------------------------------------------

bool ObsFilterData::hasVector(const std::string & grp, const std::string & var) const {
  std::map<std::string, const ioda::ObsVector *>::const_iterator jj = ovecs_.find(grp);
  if (jj == ovecs_.end()) {
    return false;
  } else {
    return jj->second->has(var);
  }
}

// -----------------------------------------------------------------------------
/*! Gets requested data from ObsFilterData
 *  \param varname is a name of a variable requested
 *  \param values on output is data from varname (undefined on input)
 *  \return data associated with varname, in std::vector<float>
 *  \warning if data are unavailable, assertions would fail and method abort
 */
void ObsFilterData::get(const Variable & varname, std::vector<float> & values) const {
  const std::string var = varname.variable();
  const std::string grp = varname.group();

  std::size_t nvals = obsdb_.nlocs();
///  VarMetaData is a special case: size(nvars) instead of (nlocs)
  if (grp == "VarMetaData")  nvals = obsdb_.nvars();

  values.resize(nvals);
///  For GeoVaLs read from GeoVaLs (should be available)
  if (grp == "GeoVaLs") {
    ASSERT(gvals_);
    gvals_->get(values, var);
///  For ObsFunction instantiate ObsFunction and calculate the result
///  TODO(AS?): cache results of function computations
  } else if (grp == "ObsFunction") {
    ioda::ObsDataVector<float> vals(obsdb_, var, grp, false);
    ObsFunction obsfunc(varname);
    obsfunc.compute(*this, vals);
    for (size_t jj = 0; jj < nvals; ++jj) {
      values[jj] = vals[var][jj];
    }
///  For HofX get from ObsVector H(x) (should be available)
  } else if (this->hasVector(grp, var)) {
    std::map<std::string, const ioda::ObsVector *>::const_iterator jv = ovecs_.find(grp);
    size_t hofxnvars = jv->second->nvars();
    size_t iv = jv->second->varnames().find(var);
    for (size_t jj = 0; jj < nvals; ++jj) {
      values[jj] = (*jv->second)[iv + (jj * hofxnvars)];
    }
///  For ObsDiag get from ObsDiagnostics
  } else if (grp == "ObsDiag") {
    ASSERT(diags_);
    diags_->get(values, var);
///  Else must be coming from ObsSpace
  } else {
    obsdb_.get_db(grp, var, values);
  }
}

// -----------------------------------------------------------------------------
/*! Gets requested integer data from ObsFilterData
 *  \param varname is a name of a variable requested
 *  \param values on output is data from varname (undefined on input)
 *  \return data associated with varname, in std::vector<int>
 *  \warning if data are unavailable, assertions would fail and method abort
 *           only ObsSpace int data are supported currently
 */
void ObsFilterData::get(const Variable & varname, std::vector<int> & values) const {
  const std::string var = varname.variable();
  const std::string grp = varname.group();

  std::size_t nvals = obsdb_.nlocs();
///  VarMetaData is a special case: size(nvars) instead of (nlocs)
  if (grp == "VarMetaData")  nvals = obsdb_.nvars();

  values.resize(nvals);
///  GeoVaLs, HofX, ObsDiag are not supportd for int data
// TODO(somebody): need something here about obs error
  if (grp == "GeoVaLs" || grp == "HofX" || grp == "ObsDiag" || grp == "ObsFunction") {
    oops::Log::error() << "ObsFilterData::get int values only supported for ObsSpace" << std::endl;
    ABORT("ObsFilterData::get int values only supported for ObsSpace");
  } else {
    obsdb_.get_db(grp, var, values);
  }
}


// -----------------------------------------------------------------------------
/*! Gets requested data at requested level from ObsFilterData
 *  \param varname is a name of a variable requested
 *         group must be either GeoVaLs or ObsDiag
 *  \param level is a level variable is requested at
 *  \param values on output is data from varname (undefined on input)
 *  \return data associated with varname, in std::vector<float>
 *  \warning if data are unavailable, assertions would fail and method abort
 */
void ObsFilterData::get(const Variable & varname, const int level,
                        std::vector<float> & values) const {
  const std::string var = varname.variable();
  const std::string grp = varname.group();

  ASSERT(grp == "GeoVaLs" || grp == "ObsDiag");
  values.resize(obsdb_.nlocs());
///  For GeoVaLs read from GeoVaLs (should be available)
  if (grp == "GeoVaLs") {
    ASSERT(gvals_);
    gvals_->get(values, var, level);
///  For ObsDiag get from ObsDiagnostics
  } else if (grp == "ObsDiag") {
    ASSERT(diags_);
    diags_->get(values, var, level);
  }
}

// -----------------------------------------------------------------------------
/*! Gets requested data from ObsFilterData into ObsDataVector
 *  \param varname is a name of a variable requested
 *  \param values on output is data from varname (should be allocated on input)
 *  \return data associated with varname, in ioda::ObsDataVector<float>
 *  \warning only works for ObsFunction;
 *           if data are unavailable, assertions would fail and method abort
 */
void ObsFilterData::get(const Variable & varname, ioda::ObsDataVector<float> & values) const {
  const std::string var = varname.variable();
  const std::string grp = varname.group();

  ASSERT(grp == "ObsFunction");

  ObsFunction obsfunc(varname);
  obsfunc.compute(*this, values);
}


// -----------------------------------------------------------------------------
/*! Returns number of levels in 3D geovals and obsdiags or
 *  one if not 3D geovals or obsdiag
 *
 */
size_t ObsFilterData::nlevs(const Variable & varname) const {
  const std::string var = varname.variable();
  const std::string grp = varname.group();
  if (grp == "GeoVaLs") {
    ASSERT(gvals_);
    return gvals_->nlevs(var);
  } else if (grp == "ObsDiag") {
    ASSERT(diags_);
    return diags_->nlevs(var);
  }
  return 1;
}

// -----------------------------------------------------------------------------
/*! Prints basic info on ObsFilterData (which data contains) */
void ObsFilterData::print(std::ostream & os) const {
  os << "Filter data: contains obs";
  if (gvals_) {
    os << ", geovals";
  }
  for (std::map<std::string, const ioda::ObsVector *>::const_iterator jj = ovecs_.begin();
       jj != ovecs_.end(); ++jj) {
    os << ", " << jj->first;
  }
  if (diags_) {
    os << ", diags";
  }
  os << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
