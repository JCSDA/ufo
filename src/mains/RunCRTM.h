/*
 * (C) Copyright 2017-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef MAINS_RUNCRTM_H_
#define MAINS_RUNCRTM_H_

#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/Observations.h"
#include "oops/base/ObsSpaces.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/interface/ObsOperator.h"
#include "oops/interface/ObsVector.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace ufo {

template <typename OBS> class RunCRTM : public oops::Application {
  typedef oops::GeoVaLs<OBS>             GeoVaLs_;
  typedef oops::ObsAuxControl<OBS>       ObsAuxCtrl_;
  typedef oops::ObsDiagnostics<OBS>      ObsDiags_;
  typedef oops::Observations<OBS>        Observations_;
  typedef oops::ObsOperator<OBS>         ObsOperator_;
  typedef oops::ObsSpaces<OBS>           ObsSpaces_;
  typedef oops::ObsVector<OBS>           ObsVector_;
  typedef oops::ObsDataVector<OBS, int>   ObsDataVectorInt_;

 public:
// -----------------------------------------------------------------------------
  explicit RunCRTM(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {}
// -----------------------------------------------------------------------------
  virtual ~RunCRTM() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig, bool validate) const {
//  Setup observation window
    const util::TimeWindow timeWindow(fullConfig.getSubConfiguration("time window"));

//  Setup observations
    ObsSpaces_ obsdb(fullConfig, this->getComm(), timeWindow);

    oops::ObsVariables diagvars;

    std::vector<eckit::LocalConfiguration> conf;
    fullConfig.get("observations", conf);

    for (std::size_t jj = 0; jj < obsdb.size(); ++jj) {
      eckit::LocalConfiguration obsopconf(conf[jj], "obs operator");
      ObsOperator_ hop(obsdb[jj], obsopconf);

      eckit::LocalConfiguration biasconf = conf[jj].getSubConfiguration("obs bias");
      const ObsAuxCtrl_ ybias(obsdb[jj], biasconf);

      oops::Variables vars = hop.requiredVars();
      oops::Variables reducedVars = ybias.requiredVars();
      vars += reducedVars;  // the reduced format is derived from the sampled format

      const eckit::LocalConfiguration gconf(conf[jj], "geovals");
      GeoVaLs_ gval(gconf, obsdb[jj], vars);
      hop.computeReducedVars(reducedVars, gval);
      ObsDataVectorInt_ qcflags(obsdb[jj], obsdb[jj].obsvariables(), std::string());
      ObsVector_ hofx(obsdb[jj]);
      ObsVector_ bias(obsdb[jj]);
      bias.zero();
      ObsDiags_ diag(obsdb[jj], hop.locations(), diagvars);

      hop.simulateObs(gval, hofx, ybias, qcflags, bias, diag);

      const double zz = hofx.rms();
      const double xx = conf[jj].getDouble("rms ref");
      const double tol = conf[jj].getDouble("tolerance");
//      BOOST_CHECK_CLOSE(xx, zz, tol);
    }

    return 0;
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const {
    return "oops::RunCRTM<" + OBS::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace ufo

#endif  // MAINS_RUNCRTM_H_
