/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/copyVars2ODV.h"

#include <string>
#include <vector>

#include "eckit/config/Configuration.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/utils/SplitVarGroup.h"

namespace ufo {


// -----------------------------------------------------------------------------

const oops::Variables collectFilterVars(const eckit::Configuration & config) {
  oops::Log::trace() << "collectFilterVars" << std::endl;

// Parse configuration
  oops::Variables allvars(config);

// where input variables
  std::vector<eckit::LocalConfiguration> masks;
  config.get("where", masks);
  for (size_t jm = 0; jm < masks.size(); ++jm) {
    allvars.push_back(masks[jm].getString("variable"));
  }

  return allvars;
}

ioda::ObsDataVector<float> copyVars2ODV(const GeoVaLs & gv,
                                        ioda::ObsDataVector<float> & odvec) {
  oops::Log::trace() << "copyVars2ODV for GeoVaLs" << std::endl;

  const oops::Variables vargrps = odvec.varnames();
  ioda::ObsDataVector<float> newODV(odvec);
  for (size_t jv = 0; jv < vargrps.size(); ++jv) {
    std::string var;
    std::string grp;
    splitVarGroup(vargrps[jv], var, grp);

    if (gv.has(var) && grp == "GeoVaLs") {
      std::vector<float> values(odvec.nlocs());

      // By default, get first level of GeoVaLs object
      gv.get(values, var);

      for (size_t jj = 0; jj < odvec.nlocs(); ++jj) {
        newODV[jv].at(jj) = values[jj];
      }
    }
  }
  return newODV;
}

ioda::ObsDataVector<float> copyVars2ODV(const ioda::ObsVector & ovec,
                                        ioda::ObsDataVector<float> & odvec,
                                        const std::string & group) {
  oops::Log::trace() << "copyVars2ODV for ObsVector" << std::endl;

  const oops::Variables vargrps = odvec.varnames();
  ioda::ObsDataVector<float> newODV(odvec);
  for (size_t jv = 0; jv < vargrps.size(); ++jv) {
    std::string var;
    std::string grp;
    splitVarGroup(vargrps[jv], var, grp);

    if (ovec.has(var) && grp == group) {
      size_t iv = ovec.varnames().find(var);
      for (size_t jj = 0; jj < odvec.nlocs(); ++jj) {
        newODV[jv].at(jj) = ovec[iv + (jj * ovec.nvars())];
      }
    }
  }
  return newODV;
}

ioda::ObsDataVector<float> copyVars2ODV(ioda::ObsSpace & obsdb,
                                        ioda::ObsDataVector<float> & odvec) {
  oops::Log::trace() << "copyVars2ODV for ObsSpace" << std::endl;

  const oops::Variables vargrps = odvec.varnames();
  ioda::ObsDataVector<float> newODV(odvec);
  for (size_t jv = 0; jv < vargrps.size(); ++jv) {
    std::string var;
    std::string grp;
    splitVarGroup(vargrps[jv], var, grp);

    if (obsdb.has(grp, var)) {
      ioda::ObsDataVector<float> vals(obsdb, var, grp);
      for (size_t jj = 0; jj < odvec.nlocs(); ++jj) {
        newODV[jv].at(jj) = vals[0].at(jj);
      }
    }
  }
  return newODV;
}

}  // namespace ufo
