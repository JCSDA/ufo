/*
 * (C) Copyright 2018-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/processWhere.h"

#include <set>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/base/Variables.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/missingValues.h"
#include "ufo/GeoVaLs.h"
#include "ufo/obsfunctions/ObsFunction.h"
#include "ufo/utils/SplitVarGroup.h"

namespace ufo {

// -----------------------------------------------------------------------------

oops::Variables preProcessWhere(const eckit::Configuration & config, const std::string & group) {
  std::vector<eckit::LocalConfiguration> masks;
  config.get("where", masks);

  std::vector<std::string> vv;
  for (size_t jm = 0; jm < masks.size(); ++jm) {
    const std::string vargrp(masks[jm].getString("variable"));
    std::string var;
    std::string grp;
    splitVarGroup(vargrp, var, grp);
    if (grp == group) vv.push_back(var);
  }

  return oops::Variables(vv);
}

// -----------------------------------------------------------------------------

bool dataExists(const std::string & varname, ioda::ObsSpace & obsdb,
                const GeoVaLs * gvals = NULL, const ioda::ObsVector * hofx = NULL) {
  std::string var, grp;
  splitVarGroup(varname, var, grp);
  if (grp == "GeoVaLs") {
    ASSERT(gvals);
    return gvals->has(var);
  } else if (grp == "ObsFunction" || grp == "HofXFunction") {
    return ObsFunctionFactory::functionExists(var);
  } else {
    return obsdb.has(grp, var);
  }
  return false;
}

// -----------------------------------------------------------------------------

std::vector<float> getData(const std::string & varname, ioda::ObsSpace & obsdb,
                           const GeoVaLs * gvals = NULL, const ioda::ObsVector * hofx = NULL) {
  std::string var, grp;
  splitVarGroup(varname, var, grp);

  std::size_t nvals = obsdb.nlocs();
  if (grp == "VarMetaData")  nvals = obsdb.nvars();

  std::vector<float> values(nvals);
  if (grp == "GeoVaLs") {
    ASSERT(gvals);
    gvals->get(values, var);
  } else if (grp == "ObsFunction") {
    ioda::ObsDataVector<float> vals(obsdb, var, grp, false);
    ObsFunction obsdiag(var);
    ioda::ObsDataVector<float> metadata(obsdb, obsdiag.requiredMetaData(), "MetaData");
    ioda::ObsDataVector<float> obs(obsdb, obsdiag.requiredObsData(), "ObsValue");
    obsdiag.compute(metadata, obs, *gvals, vals);
    for (size_t jj = 0; jj < nvals; ++jj) {
      values[jj] = vals[var][jj];
    }
  } else if (grp == "HofXFunction") {
    ASSERT(hofx);
    ioda::ObsDataVector<float> vals(obsdb, var, grp, false);
    ObsFunction obsdiag(var);
    ioda::ObsDataVector<float> metadata(obsdb, obsdiag.requiredMetaData(), "MetaData");
    const oops::Variables requiredObs = obsdiag.requiredObsData();
    ioda::ObsDataVector<float> obs(obsdb, obsdiag.requiredObsData(), "HofX", false);
    const size_t hofxnvars = hofx->nvars();
    for (size_t jv = 0; jv < requiredObs.size(); ++jv) {
      ASSERT(hofx->has(requiredObs[jv]));
      size_t iv = hofx->varnames().find(requiredObs[jv]);
      for (size_t jj = 0; jj < obs.nlocs(); ++jj) {
        obs[jv].at(jj) = (*hofx)[iv + (jj * hofxnvars)];
      }
    }
    obsdiag.compute(metadata, obs, *gvals, vals);
    for (size_t jj = 0; jj < nvals; ++jj) {
      values[jj] = vals[var][jj];
    }
  } else if (grp == "HofX") {
    ASSERT(hofx);
    ASSERT(hofx->has(var));
    size_t hofxnvars = hofx->nvars();
    size_t iv = hofx->varnames().find(var);
    for (size_t jj = 0; jj < nvals; ++jj) {
      values[jj] = (*hofx)[iv + (jj * hofxnvars)];
    }
  } else {
    obsdb.get_db(grp, var, nvals, values.data());
  }
  return values;
}

// -----------------------------------------------------------------------------

void processWhereMinMax(const std::vector<float> & data,
            const float & vmin, const float & vmax, std::vector<bool> & mask) {
  const float missing = util::missingValue(missing);
  const size_t n = data.size();

  if (vmin != missing || vmax != missing) {
    for (size_t jj = 0; jj < n; ++jj) {
      if (vmin != missing && data[jj] < vmin) mask[jj] = false;
      if (vmax != missing && data[jj] > vmax) mask[jj] = false;
    }
  }
}

// -----------------------------------------------------------------------------

void processWhereIsDefined(const std::vector<float> & data,
                           std::vector<bool> & mask) {
  const float missing = util::missingValue(missing);
  const size_t n = data.size();
  for (size_t jj = 0; jj < n; ++jj) {
    if (data[jj] == missing) mask[jj] = false;
  }
}

// -----------------------------------------------------------------------------

void processWhereIsNotDefined(const std::vector<float> & data,
                              std::vector<bool> & mask) {
  const float missing = util::missingValue(missing);
  const size_t n = data.size();
  for (size_t jj = 0; jj < n; ++jj) {
    if (data[jj] != missing) mask[jj] = false;
  }
}

// -----------------------------------------------------------------------------

void processWhereIsIn(const std::vector<int> & data,
                      const std::set<int> & whitelist,
                      std::vector<bool> & mask) {
  const size_t n = data.size();
  for (size_t jj = 0; jj < n; ++jj) {
    if (!oops::contains(whitelist, data[jj])) mask[jj] = false;
  }
}

// -----------------------------------------------------------------------------

void processWhereIsNotIn(const std::vector<int> & data,
                         const std::set<int> & blacklist,
                         std::vector<bool> & mask) {
  const size_t n = data.size();
  for (size_t jj = 0; jj < n; ++jj) {
    if (oops::contains(blacklist, data[jj])) mask[jj] = false;
  }
}

// -----------------------------------------------------------------------------

std::vector<bool> processWhere(const eckit::Configuration & config,
        ioda::ObsSpace & obsdb, const GeoVaLs * gvals,
        const ioda::ObsVector * hofx) {
  const float missing = util::missingValue(missing);
  const size_t nlocs = obsdb.nlocs();

// Everywhere by default if no mask
  std::vector<bool> where(nlocs, true);

  std::vector<eckit::LocalConfiguration> masks;
  config.get("where", masks);

  for (size_t jm = 0; jm < masks.size(); ++jm) {
//  Get variable@group
    const std::string varname(masks[jm].getString("variable"));
    std::string var, grp;
    splitVarGroup(varname, var, grp);
    if (grp != "VarMetaData") {
//    Get data
      std::vector<float> data = getData(varname, obsdb, gvals, hofx);
//    Process masks on float values
      const float vmin = masks[jm].getFloat("minvalue", missing);
      const float vmax = masks[jm].getFloat("maxvalue", missing);

//    Apply mask min/max
      if (vmin != missing || vmax != missing) {
        processWhereMinMax(data, vmin, vmax, where);
      }

//    Apply mask is_defined
      if (masks[jm].has("is_defined")) {
        if (dataExists(varname, obsdb, gvals, hofx)) {
          processWhereIsDefined(data, where);
        } else {
          std::fill(where.begin(), where.end(), false);
        }
      }

//    Apply mask is_not_defined
      if (masks[jm].has("is_not_defined")) {
        processWhereIsNotDefined(data, where);
      }
    }
  }
//  Print diagnostics for debug
  int ii = 0;
  for (size_t jj = 0; jj < nlocs; ++jj) {
    if (where[jj] == false) ++ii;
  }

  oops::Log::debug() << "processWhere: selected " << ii << " obs." << std::endl;
  return where;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
