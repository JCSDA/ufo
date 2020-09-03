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
  : obsdb_(obsdb), gvals_(NULL), ovecs_(), diags_(NULL), dvecsf_(), dvecsi_() {
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
/*! Associates ObsDataVector with this ObsFilterData */
void ObsFilterData::associate(const ioda::ObsDataVector<float> & data, const std::string & name) {
  dvecsf_[name] = &data;
}

// -----------------------------------------------------------------------------
/*! Associates ObsDataVector with this ObsFilterData */
void ObsFilterData::associate(const ioda::ObsDataVector<int> & data, const std::string & name) {
  dvecsi_[name] = &data;
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
  } else if (grp == "ObsDiag" || grp == "ObsBiasTerm") {
    return (diags_ && diags_->has(var));
  } else {
    return this->hasVector(grp, var) || this->hasDataVector(grp, var) || obsdb_.has(grp, var);
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

bool ObsFilterData::hasDataVector(const std::string & grp, const std::string & var) const {
  std::map<std::string, const ioda::ObsDataVector<float> *>::const_iterator jj = dvecsf_.find(grp);
  if (jj == dvecsf_.end()) {
    return false;
  } else {
    return jj->second->has(var);
  }
}

// -----------------------------------------------------------------------------

bool ObsFilterData::hasDataVectorInt(const std::string & grp, const std::string & var) const {
  std::map<std::string, const ioda::ObsDataVector<int> *>::const_iterator jj = dvecsi_.find(grp);
  if (jj == dvecsi_.end()) {
    return false;
  } else {
    return jj->second->has(var);
  }
}

// -----------------------------------------------------------------------------
/*! Gets requested data from ObsFilterData
 *  \param[in] varname is a name of a variable requested
 *  \param[out] values on output is data from varname (undefined on input)
 *  \warning if data are unavailable, assertions would fail and method abort
 */
void ObsFilterData::get(const Variable & varname, std::vector<float> & values) const {
  const std::string var = varname.variable(0);
  const std::string grp = varname.group();

  if (grp == "VarMetaData") {
    values.resize(obsdb_.nvars());
    obsdb_.get_db(grp, var, values);
  } else {
    ioda::ObsDataVector<float> vec(obsdb_, varname.toOopsVariables(), grp, false);
    this->get(varname, vec);
    values.resize(obsdb_.nlocs());
    for (size_t jj = 0; jj < obsdb_.nlocs(); ++jj) {
      values[jj] = vec[var][jj];
    }
  }
}


// -----------------------------------------------------------------------------
/*! Gets requested data from ObsFilterData
 *  \param[in] varname is a name of a variable requested
 *  \param[out] values on output is data from varname (undefined on input)
 *  \warning if data are unavailable, assertions would fail and method abort
 */
void ObsFilterData::get(const Variable & varname, std::vector<std::string> & values) const {
  const std::string var = varname.variable();
  const std::string grp = varname.group();

  if (grp == "GeoVaLs" || grp == "HofX" || grp == "ObsDiag" ||
    grp == "ObsBiasTerm" || grp == "ObsFunction") {
      oops::Log::error() << "ObsFilterData::get std::string and int values only "
                         << "supported for ObsSpace"
                         << std::endl;
    ABORT("ObsFilterData::get std::string, int and util::DateTime values only supported for "
          "ObsSpace");
  } else {
    values.resize(obsdb_.nlocs());
    obsdb_.get_db(grp, var, values);
  }
}


// -----------------------------------------------------------------------------
/*! Gets requested data from ObsFilterData
 *  \param[in] varname is a name of a variable requested
 *  \param[out] values on output is data from varname (undefined on input)
 *  \warning if data are unavailable, assertions would fail and method abort
 */
void ObsFilterData::get(const Variable & varname, std::vector<util::DateTime> & values) const {
  const std::string var = varname.variable();
  const std::string grp = varname.group();

  if (grp == "GeoVaLs" || grp == "HofX" || grp == "ObsDiag" ||
    grp == "ObsBiasTerm" || grp == "ObsFunction") {
      oops::Log::error() << "ObsFilterData::get std::string and int values only "
                         << "supported for ObsSpace"
                         << std::endl;
    ABORT("ObsFilterData::get std::string, int and util::DateTime values only supported for "
          "ObsSpace");
  } else {
    values.resize(obsdb_.nlocs());
    obsdb_.get_db(grp, var, values);
  }
}


// -----------------------------------------------------------------------------
/*! Gets requested integer data from ObsFilterData
 *  \param[in] varname is a name of a variable requested
 *  \param[out] values on output is data from varname (undefined on input)
 *  \warning if data are unavailable, assertions would fail and method abort
 *           only ObsSpace int data are supported currently
 */
void ObsFilterData::get(const Variable & varname, std::vector<int> & values) const {
  const std::string var = varname.variable();
  const std::string grp = varname.group();

  if (grp == "VarMetaData") {
    values.resize(obsdb_.nvars());
    obsdb_.get_db(grp, var, values);
  } else {
    if (grp == "GeoVaLs" || grp == "HofX" || grp == "ObsDiag" ||
        grp == "ObsBiasTerm" || grp == "ObsFunction") {
      oops::Log::error() << "ObsFilterData::get std::string and int values only "
                         << "supported for ObsSpace"
                         << std::endl;
      ABORT("ObsFilterData::get std::string and int values only supported for ObsSpace");
    } else {
      ioda::ObsDataVector<int> vec(obsdb_, varname.toOopsVariables(), grp, false);
      this->get(varname, vec);
      values.resize(obsdb_.nlocs());
      for (size_t jj = 0; jj < obsdb_.nlocs(); ++jj) {
        values[jj] = vec[var][jj];
      }
    }
  }
}


// -----------------------------------------------------------------------------
/*! Gets requested data at requested level from ObsFilterData
 *  \param[in] varname is a name of a variable requested
 *         group must be either GeoVaLs or ObsDiag
 *  \param[in] level is a level variable is requested at
 *  \param[out] values on output is data from varname (undefined on input)
 *  \warning if data are unavailable, assertions would fail and method abort
 */
void ObsFilterData::get(const Variable & varname, const int level,
                        std::vector<float> & values) const {
  const std::string var = varname.variable();
  const std::string grp = varname.group();

  ASSERT(grp == "GeoVaLs" || grp == "ObsDiag" || grp == "ObsBiasTerm");
  values.resize(obsdb_.nlocs());
///  For GeoVaLs read from GeoVaLs (should be available)
  if (grp == "GeoVaLs") {
    ASSERT(gvals_);
    gvals_->get(values, var, level);
///  For ObsDiag get from ObsDiagnostics
  } else if (grp == "ObsDiag" || grp == "ObsBiasTerm") {
    ASSERT(diags_);
    diags_->get(values, var, level);
  }
}

// -----------------------------------------------------------------------------
/*! Gets requested data from ObsFilterData into ObsDataVector
 *  \param[in] varname is a name of a variable requested
 *  \param[out] values on output is data from varname (should be allocated on input)
 *  \warning only works for ObsFunction;
 *           if data are unavailable, assertions would fail and method abort
 */
void ObsFilterData::get(const Variable & varname, ioda::ObsDataVector<float> & values) const {
  const std::string var = varname.variable(0);
  const std::string grp = varname.group();
  /// For GeoVaLs read single variable and save in the relevant field
  if (grp == "GeoVaLs") {
    ASSERT(gvals_);
    std::vector<float> vec(obsdb_.nlocs());
    gvals_->get(vec, var);
    values[var] = vec;
  /// For Function call compute
  } else if (grp == "ObsFunction") {
    ObsFunction obsfunc(varname);
    obsfunc.compute(*this, values);
  ///  For HofX get from ObsVector H(x) (should be available)
  } else if (this->hasVector(grp, var)) {
    std::map<std::string, const ioda::ObsVector *>::const_iterator jv = ovecs_.find(grp);
    size_t hofxnvars = jv->second->nvars();
    for (size_t ivar = 0; ivar < varname.size(); ++ivar) {
      const std::string currvar = varname.variable(ivar);
      size_t iv = jv->second->varnames().find(currvar);
      for (size_t jj = 0; jj < obsdb_.nlocs(); ++jj) {
        values[currvar][jj] = (*jv->second)[iv + (jj * hofxnvars)];
      }
    }
///  For ObsDiag or ObsBiasTerm,  get it from ObsDiagnostics
  } else if (grp == "ObsDiag" || grp == "ObsBiasTerm") {
    ASSERT(diags_);
    std::vector<float> vec(obsdb_.nlocs());
    diags_->get(vec, var);
    values[var] = vec;
///  For ObsDataVector
  } else if (this->hasDataVector(grp, var)) {
    std::map<std::string, const ioda::ObsDataVector<float> *>::const_iterator
                          jv = dvecsf_.find(grp);
    values = *jv->second;
  } else {
    values.read(grp);
  }
}

// -----------------------------------------------------------------------------
/*! Gets requested data from ObsFilterData into ObsDataVector
 *  \param[in] varname is a name of a variable requested
 *  \param[out] values on output is data from varname (should be allocated on input)
 *  \return data associated with varname, in ioda::ObsDataVector<int>
*/
void ObsFilterData::get(const Variable & varname, ioda::ObsDataVector<int> & values) const {
  const std::string var = varname.variable(0);
  const std::string grp = varname.group();
  if (this->hasDataVectorInt(grp, var)) {
    std::map<std::string, const ioda::ObsDataVector<int> *>::const_iterator jv = dvecsi_.find(grp);
    values = *jv->second;
  } else {
    values.read(grp);
  }
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
  } else if (grp == "ObsDiag" || grp == "ObsBiasTerm") {
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
  for (std::map<std::string, const ioda::ObsDataVector<float> *>::const_iterator
     jj = dvecsf_.begin(); jj != dvecsf_.end(); ++jj) {
    os << ", " << jj->first;
  }
  for (std::map<std::string, const ioda::ObsDataVector<int> *>::const_iterator
     jj = dvecsi_.begin(); jj != dvecsi_.end(); ++jj) {
    os << ", " << jj->first;
  }
  if (diags_) {
    os << ", diags";
  }
  os << std::endl;
}

// -----------------------------------------------------------------------------
/*! \brief dtype of the provided variable
 *  \details This method returns the data type of the variable stored in the
 *           ObsFilterData.
 *  \param varname is a name of a variable requested
 *  \return data type (ioda::ObsDtype) associated with varname
 */
ioda::ObsDtype ObsFilterData::dtype(const Variable & varname) const {
  const std::string var = varname.variable();
  const std::string grp = varname.group();
  // Default to float
  ioda::ObsDtype res = ioda::ObsDtype::Float;
  if (obsdb_.has(grp, var)) {
    res = obsdb_.dtype(grp, var);
  } else if (!this->has(varname)) {
    oops::Log::error() << "ObsFilterData::dtype unable to find provided variable."
                       << std::endl;
    ABORT("ObsFilterData::dtype unable to find provided variable.");
  }
  return res;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
