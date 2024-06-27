/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/ObsFilterData.h"

#include <string>
#include <vector>

#include "eckit/utils/StringTools.h"
#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/obsfunctions/ObsFunction.h"
#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variable.h"
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
/*! Associates H(x)-like ObsVector with this ObsFilterData */
void ObsFilterData::associate(const ioda::ObsVector & data, const std::string & name) {
  ovecs_[name] = &data;
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
    return (gvals_ && gvals_->has(oops::Variable(var)));
  } else if (grp == ObsFunctionTraits<float>::groupName) {
    return ObsFunctionFactory<float>::functionExists(var);
  } else if (grp == ObsFunctionTraits<int>::groupName) {
    return ObsFunctionFactory<int>::functionExists(var);
  } else if (grp == ObsFunctionTraits<std::string>::groupName) {
    return ObsFunctionFactory<std::string>::functionExists(var);
  } else if (grp == ObsFunctionTraits<util::DateTime>::groupName) {
    return ObsFunctionFactory<util::DateTime>::functionExists(var);
  } else if (grp == "ObsDiag" || grp == "ObsBiasTerm") {
    return (diags_ && diags_->has(var));
  } else {
    return this->hasVector(grp, var) ||
           this->hasDataVector(grp, var) ||
           this->hasDataVectorInt(grp, var) ||
           obsdb_.has(grp, var);
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
// Overloads of get() taking a std::vector.

// -----------------------------------------------------------------------------
void ObsFilterData::get(const Variable & varname, std::vector<float> & values,
                        bool skipDerived) const {
  getVector(varname, values, skipDerived);
}

// -----------------------------------------------------------------------------
void ObsFilterData::get(const Variable & varname, std::vector<int> & values,
                        bool skipDerived) const {
  getVector(varname, values, skipDerived);
}

// -----------------------------------------------------------------------------
void ObsFilterData::get(const Variable & varname, std::vector<std::string> & values,
                        bool skipDerived) const {
  getVector(varname, values, skipDerived);
}

// -----------------------------------------------------------------------------
void ObsFilterData::get(const Variable & varname, std::vector<util::DateTime> & values,
                        bool skipDerived) const {
  getVector(varname, values, skipDerived);
}

// -----------------------------------------------------------------------------
void ObsFilterData::get(const Variable & varname, std::vector<DiagnosticFlag> & values,
                        bool skipDerived) const {
  getVector(varname, values, skipDerived);
}

// -----------------------------------------------------------------------------
template <typename T>
void ObsFilterData::getVector(const Variable & varname, std::vector<T> & values,
                              bool skipDerived) const {
  const std::string var = varname.variable();
  const std::string grp = varname.group();
  if (obsdb_.has(grp, var, skipDerived)) {
    obsdb_.get_db(grp, var, values, {}, skipDerived);
  } else {
    ioda::ObsDataVector<T> vec(obsdb_, varname.toOopsObsVariables(), grp, false);
    this->get(varname, vec, skipDerived);
    values = vec[var];
  }
}

// -----------------------------------------------------------------------------
// Overload of get() taking an std::vector<float> and a level index.
// -----------------------------------------------------------------------------
void ObsFilterData::get(const Variable & varname, const int level,
                        std::vector<float> & values) const {
  const std::string var = varname.variable();
  const std::string grp = varname.group();

  ASSERT(grp == "GeoVaLs" || grp == "ObsDiag" || grp == "ObsBiasTerm");
  values.resize(obsdb_.nlocs());
///  For GeoVaLs read from GeoVaLs (should be available)
  if (grp == "GeoVaLs") {
    ASSERT(gvals_);
    gvals_->getAtLevel(values, oops::Variable(var), level);
///  For ObsDiag get from ObsDiagnostics
  } else if (grp == "ObsDiag" || grp == "ObsBiasTerm") {
    ASSERT(diags_);
    diags_->get(values, var, level);
  }
}

// -----------------------------------------------------------------------------
// Overload of get() taking an std::vector<double> and a level index.
// -----------------------------------------------------------------------------
void ObsFilterData::get(const Variable & varname, const int level,
                        std::vector<double> & values) const {
  const std::string var = varname.variable();
  const std::string grp = varname.group();

  ASSERT(grp == "GeoVaLs" || grp == "ObsDiag" || grp == "ObsBiasTerm");
  values.resize(obsdb_.nlocs());
///  For GeoVaLs read from GeoVaLs (should be available)
  if (grp == "GeoVaLs") {
    ASSERT(gvals_);
    gvals_->getAtLevel(values, oops::Variable(var), level);
///  For ObsDiag get from ObsDiagnostics
  } else if (grp == "ObsDiag" || grp == "ObsBiasTerm") {
    ASSERT(diags_);
    diags_->get(values, var, level);
  }
}

// -----------------------------------------------------------------------------
// Overloads of get() taking an ioda::ObsDataVector.

// -----------------------------------------------------------------------------
void ObsFilterData::get(const Variable & varname, ioda::ObsDataVector<float> & values,
                        bool skipDerived) const {
  const std::string var = varname.variable(0);
  const std::string grp = varname.group();
  /// For GeoVaLs read single variable and save in the relevant field
  if (grp == "GeoVaLs") {
    ASSERT(gvals_);
    std::vector<float> vec(obsdb_.nlocs());
    gvals_->get(vec, oops::Variable(var));
    values[var] = vec;
  /// For Function call compute
  } else if (grp == ObsFunctionTraits<float>::groupName) {
    ObsFunction<float> obsfunc(varname);
    obsfunc.compute(*this, values);
  ///  For HofX get from ObsVector H(x) (should be available)
  } else if (this->hasVector(grp, var)) {
    std::map<std::string, const ioda::ObsVector *>::const_iterator jv = ovecs_.find(grp);
    // copy H(x) ObsVector to 'values', converting missing vals from double to float:
    values.assignToExistingVariables(*jv->second);
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
    for (size_t ivar = 0; ivar < varname.size(); ++ivar) {
      const std::string currvar = varname.variable(ivar);
      values[currvar] = (*jv->second)[currvar];
    }
  /// Produce a clear error message for incompatible ObsFunction types
  } else if (eckit::StringTools::endsWith(grp, "ObsFunction")) {
    throw eckit::BadParameter("ObsFilterData::get(): " + varname.fullName() +
                              " is not a function producing values of type " +
                              ObsFunctionTraits<float>::valueTypeName, Here());
  } else {
    values.read(grp, true, skipDerived);
  }
}

// -----------------------------------------------------------------------------
void ObsFilterData::get(const Variable & varname, ioda::ObsDataVector<int> & values,
                        bool skipDerived) const {
  const std::string var = varname.variable(0);
  const std::string grp = varname.group();
  /// For Function call compute
  if (grp == ObsFunctionTraits<int>::groupName) {
    ObsFunction<int> obsfunc(varname);
    obsfunc.compute(*this, values);
  /// For ObsDataVector
  } else if (this->hasDataVectorInt(grp, var)) {
    std::map<std::string, const ioda::ObsDataVector<int> *>::const_iterator jv = dvecsi_.find(grp);
    for (size_t ivar = 0; ivar < varname.size(); ++ivar) {
      const std::string currvar = varname.variable(ivar);
      values[currvar] = (*jv->second)[currvar];
    }
  /// Produce a clear error message for incompatible ObsFunction types
  } else if (eckit::StringTools::endsWith(grp, "ObsFunction")) {
    throw eckit::BadParameter("ObsFilterData::get(): " + varname.fullName() +
                              " is not a function producing values of type " +
                              ObsFunctionTraits<int>::valueTypeName, Here());
  /// Produce a clear error message for other incompatible groups
  } else if (grp == "GeoVaLs" || grp == "HofX" || grp == "ObsDiag" || grp == "ObsBiasTerm") {
    throw eckit::BadParameter("ObsFilterData::get(): variables from the group " + grp +
                              " are of type float", Here());
  } else {
    values.read(grp, true, skipDerived);
  }
}

// -----------------------------------------------------------------------------
void ObsFilterData::get(const Variable & varname, ioda::ObsDataVector<std::string> & values,
                        bool skipDerived) const {
  getNonNumeric(varname, values, skipDerived);
}

// -----------------------------------------------------------------------------
void ObsFilterData::get(const Variable & varname,
                        ioda::ObsDataVector<util::DateTime> & values,
                        bool skipDerived) const {
  getNonNumeric(varname, values, skipDerived);
}

// -----------------------------------------------------------------------------
void ObsFilterData::get(const Variable & varname,
                        ioda::ObsDataVector<DiagnosticFlag> & values,
                        bool skipDerived) const {
  const std::string &grp = varname.group();
  // There are no ObsFunctions producing flags yet
  if (eckit::StringTools::endsWith(grp, "ObsFunction")) {
    throw eckit::BadParameter("ObsFilterData::get(): " + varname.fullName() +
                              " does not produce values of type DiagnosticFlag", Here());
  // Certain other groups can't contain diagnostic flags either
  } else if (grp == "GeoVaLs" || grp == "HofX" || grp == "ObsDiag" || grp == "ObsBiasTerm") {
    throw eckit::BadParameter("ObsFilterData::get(): variables from the group " + grp +
                              " are not of type DiagnosticFlag", Here());
  } else {
    // Try to retrieve the requested diagnostic flag from the ObsSpace. read() will throw an
    // exception if that flag is not found.
    values.read(grp, true, skipDerived);
  }
}

// -----------------------------------------------------------------------------
template <typename T>
void ObsFilterData::getNonNumeric(const Variable & varname, ioda::ObsDataVector<T> & values,
                                  bool skipDerived) const {
  const std::string &grp = varname.group();
  /// For Function call compute
  if (grp == ObsFunctionTraits<T>::groupName) {
    ObsFunction<T> obsfunc(varname);
    obsfunc.compute(*this, values);
  } else if (eckit::StringTools::endsWith(grp, "ObsFunction")) {
    throw eckit::BadParameter("ObsFilterData::get(): " + varname.fullName() +
                              " is not a function producing values of type " +
                              ObsFunctionTraits<T>::valueTypeName, Here());
  /// Produce a clear error message for other incompatible groups
  } else if (grp == "GeoVaLs" || grp == "HofX" || grp == "ObsDiag" || grp == "ObsBiasTerm") {
    throw eckit::BadParameter("ObsFilterData::get(): variables from the group " + grp +
                              " are of type float", Here());
  } else {
    values.read(grp, true, skipDerived);
  }
}

// -----------------------------------------------------------------------------
// End of overloads of get().

// -----------------------------------------------------------------------------
size_t ObsFilterData::nlevs(const Variable & varname) const {
  const std::string var = varname.variable();
  const std::string grp = varname.group();
  if (grp == "GeoVaLs") {
    ASSERT(gvals_);
    return gvals_->nlevs(oops::Variable(var));
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

  // Print quantities in more detail
  os << "Contents of ObsValue" << std::endl;
  for (const std::string & varname : obsdb_.obsvariables().variables()) {
    os << "- " << varname << std::endl;
  }
  os << std::endl;
  if (gvals_) {
    os << "Contents of GeoVaLs" << std::endl;
    for (const auto & var : gvals_->getVars()) {
      os << "- " << var.name() << std::endl;
    }
    os << std::endl;
  }
  for (const auto & ovec : ovecs_) {
    os << "Contents of " << ovec.first << std::endl;
    for (const std::string & varname : ovec.second->varnames().variables())
      os << "- " << varname << std::endl;
    os << std::endl;
  }
  for (const auto & dvecf : dvecsf_) {
    os << "Contents of " << dvecf.first << std::endl;
    for (const std::string & varname : dvecf.second->varnames().variables())
      os << "- " << varname << std::endl;
    os << std::endl;
  }
  for (const auto & dveci : dvecsi_) {
    os << "Contents of " << dveci.first << std::endl;
    for (const std::string & varname : dveci.second->varnames().variables())
      os << "- " << varname << std::endl;
    os << std::endl;
  }
  // todo: add accessor method to contents of ObsDiagnostics
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
  if (obsdb_.empty()) {
    res = ioda::ObsDtype::Empty;
  } else if (obsdb_.has(grp, var)) {
    res = obsdb_.dtype(grp, var);
  } else if (hasDataVectorInt(varname.group(), varname.variable())) {
    res = ioda::ObsDtype::Integer;
  } else if (grp == ObsFunctionTraits<int>::groupName &&
             ObsFunctionFactory<int>::functionExists(var)) {
    res = ioda::ObsDtype::Integer;
  } else if (grp == ObsFunctionTraits<std::string>::groupName &&
             ObsFunctionFactory<std::string>::functionExists(var)) {
    res = ioda::ObsDtype::String;
  } else if (grp == ObsFunctionTraits<util::DateTime>::groupName &&
             ObsFunctionFactory<util::DateTime>::functionExists(var)) {
    res = ioda::ObsDtype::DateTime;
  } else if (!this->has(varname)) {
    oops::Log::error() << "ObsFilterData::dtype unable to find provided variable: "
                       << varname.fullName()
                       << std::endl;
    ABORT("ObsFilterData::dtype unable to find provided variable.");
  }
  return res;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
