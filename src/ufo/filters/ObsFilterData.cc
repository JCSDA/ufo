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
#include "ufo/utils/SplitVarGroup.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsFilterData::ObsFilterData(ioda::ObsSpace & obsdb)
  : obsdb_(obsdb), gvals_(), hofx_() {
  oops::Log::trace() << "ObsFilterData created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsFilterData::~ObsFilterData() {
  oops::Log::trace() << "ObsFilterData destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsFilterData::associate(const GeoVaLs & gvals) {
  gvals_ = &gvals;
}

// -----------------------------------------------------------------------------

void ObsFilterData::associate(const ioda::ObsVector & hofx) {
  hofx_ = &hofx;
}

// -----------------------------------------------------------------------------

size_t ObsFilterData::nlocs() const {
  return obsdb_.nlocs();
}

// -----------------------------------------------------------------------------

bool ObsFilterData::has(const std::string & varname) const {
  std::string var, grp;
  splitVarGroup(varname, var, grp);
  if (grp == "GeoVaLs") {
    ASSERT(gvals_);
    return gvals_->has(var);
  } else if (grp == "ObsFunction" || grp == "HofXFunction") {
    return ObsFunctionFactory::functionExists(var);
  } else {
    return obsdb_.has(grp, var);
  }
  return false;
}

// -----------------------------------------------------------------------------

std::vector<float> ObsFilterData::get(const std::string & varname) const {
  std::string var, grp;
  splitVarGroup(varname, var, grp);

  std::size_t nvals = obsdb_.nlocs();
  if (grp == "VarMetaData")  nvals = obsdb_.nvars();

  std::vector<float> values(nvals);
  if (grp == "GeoVaLs") {
    ASSERT(gvals_);
    gvals_->get(values, var);
  } else if (grp == "ObsFunction" || grp == "HofXFunction") {
    ioda::ObsDataVector<float> vals(obsdb_, var, grp, false);
    ObsFunction obsdiag(var);
    obsdiag.compute(*this, vals);
    for (size_t jj = 0; jj < nvals; ++jj) {
      values[jj] = vals[var][jj];
    }
  } else if (grp == "HofX") {
    ASSERT(hofx_);
    ASSERT(hofx_->has(var));
    size_t hofxnvars = hofx_->nvars();
    size_t iv = hofx_->varnames().find(var);
    for (size_t jj = 0; jj < nvals; ++jj) {
      values[jj] = (*hofx_)[iv + (jj * hofxnvars)];
    }
  } else {
    obsdb_.get_db(grp, var, nvals, values.data());
  }
  return values;
}

// -----------------------------------------------------------------------------

void ObsFilterData::print(std::ostream & os) const {
  os << "Filter data: contains obs";
  if (gvals_) {
    os << ", geovals";
  }
  if (hofx_) {
    os << ", hofx";
  }
  os << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
