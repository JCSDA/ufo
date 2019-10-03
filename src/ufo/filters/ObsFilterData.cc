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
#include "oops/util/IntSetParser.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/obsfunctions/ObsFunction.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/utils/SplitVarGroup.h"

namespace ufo {

// -----------------------------------------------------------------------------
ObsFilterData::ObsFilterData(ioda::ObsSpace & obsdb)
  : obsdb_(obsdb), gvals_(NULL), hofx_(NULL), diags_(NULL) {
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
void ObsFilterData::associate(const ioda::ObsVector & hofx) {
  hofx_ = &hofx;
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
 *  \param varname is a name of a variable requested (has to be formatted as
 *         name@group
 *  \return true if the variable is available for access from ObsFilterData,
 *          false otherwise
 */
bool ObsFilterData::has(const std::string & varname) const {
  std::string var, grp;
  splitVarGroup(varname, var, grp);
  if (grp == "GeoVaLs") {
    return (gvals_ && gvals_->has(var));
  } else if (grp == "ObsFunction") {
    return ObsFunctionFactory::functionExists(var);
  } else if (grp == "HofX") {
    return (hofx_ && hofx_->has(var));
  } else if (grp == "ObsDiag") {
    return (diags_ && diags_->has(var));
  } else {
    return obsdb_.has(grp, var);
  }
  return false;
}

// -----------------------------------------------------------------------------
/*! Gets requested data from ObsFilterData
 *  \param varname is a name of a variable requested (has to be formatted as
 *         name@group
 *  \param values on output is data from varname (undefined on input)
 *  \return data associated with varname, in std::vector<float>
 *  \warning if data are unavailable, assertions would fail and method abort
 */
void ObsFilterData::get(const std::string & varname, std::vector<float> & values) const {
  std::string var, grp;
  splitVarGroup(varname, var, grp);

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
    ObsFunction obsfunc(var);
    obsfunc.compute(*this, vals);
    for (size_t jj = 0; jj < nvals; ++jj) {
      values[jj] = vals[var][jj];
    }
///  For HofX get from ObsVector H(x) (should be available)
  } else if (grp == "HofX") {
    ASSERT(hofx_);
    ASSERT(hofx_->has(var));
    size_t hofxnvars = hofx_->nvars();
    size_t iv = hofx_->varnames().find(var);
    for (size_t jj = 0; jj < nvals; ++jj) {
      values[jj] = (*hofx_)[iv + (jj * hofxnvars)];
    }
///  For ObsDiag get from ObsDiagnostics
  } else if (grp == "ObsDiag") {
    ASSERT(diags_);
    diags_->get(values, var);
///  Else must be coming from ObsSpace
  } else {
    obsdb_.get_db(grp, var, nvals, values.data());
  }
}

// -----------------------------------------------------------------------------
/*! Gets requested integer data from ObsFilterData
 *  \param varname is a name of a variable requested (has to be formatted as
 *         name@group
 *  \param values on output is data from varname (undefined on input)
 *  \return data associated with varname, in std::vector<int>
 *  \warning if data are unavailable, assertions would fail and method abort
 *           only ObsSpace int data are supported currently
 */
void ObsFilterData::get(const std::string & varname, std::vector<int> & values) const {
  std::string var, grp;
  splitVarGroup(varname, var, grp);

  std::size_t nvals = obsdb_.nlocs();
///  VarMetaData is a special case: size(nvars) instead of (nlocs)
  if (grp == "VarMetaData")  nvals = obsdb_.nvars();

  values.resize(nvals);
///  GeoVaLs, HofX, ObsDiag are not supportd for int data
  if (grp == "GeoVaLs" || grp == "HofX" || grp == "ObsDiag" || grp == "ObsFunction") {
    oops::Log::error() << "ObsFilterData::get int values only supported for ObsSpace" << std::endl;
    ABORT("ObsFilterData::get int values only supported for ObsSpace");
  } else {
    obsdb_.get_db(grp, var, nvals, values.data());
  }
}


// -----------------------------------------------------------------------------
/*! Gets requested data at requested level from ObsFilterData
 *  \param varname is a name of a variable requested (has to be formatted as
 *         name@group, group must be either GeoVaLs or ObsDiag
 *  \param level is a level variable is requested at
 *  \param values on output is data from varname (undefined on input)
 *  \return data associated with varname, in std::vector<float>
 *  \warning if data are unavailable, assertions would fail and method abort
 */
void ObsFilterData::get(const std::string & varname, const int level,
                        std::vector<float> & values) const {
  std::string var, grp;
  splitVarGroup(varname, var, grp);

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
 *  \param varname is a name of a variable requested (has to be formatted as
 *         name@group
 *  \param values on output is data from varname (should be allocated on input)
 *  \return data associated with varname, in ioda::ObsDataVector<float>
 *  \warning only works for ObsFunction;
 *           if data are unavailable, assertions would fail and method abort
 */
void ObsFilterData::get(const std::string & varname, ioda::ObsDataVector<float> & values) const {
  std::string var, grp;
  splitVarGroup(varname, var, grp);

  ASSERT(grp == "ObsFunction");

  ObsFunction obsfunc(var);
  obsfunc.compute(*this, values);
}


// -----------------------------------------------------------------------------
/*! Returns number of levels in 3D geovals and obsdiags or
 *  one if not 3D geovals or obsdiag
 *
 */
size_t ObsFilterData::nlevs(const std::string & varname) const {
  std::string var, grp;
  splitVarGroup(varname, var, grp);
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
  if (hofx_) {
    os << ", hofx";
  }
  if (diags_) {
    os << ", diags";
  }
  os << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
