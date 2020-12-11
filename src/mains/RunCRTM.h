/*
 * (C) Copyright 2017 UCAR
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
#include "oops/interface/ObsOperator.h"
#include "oops/interface/ObsVector.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace ufo {

template <typename MODEL> class RunCRTM : public oops::Application {
  typedef oops::GeoVaLs<MODEL>           GeoVaLs_;
  typedef oops::ObsAuxControl<MODEL>     ObsAuxCtrl_;
  typedef oops::ObsDiagnostics<MODEL>    ObsDiags_;
  typedef oops::Observations<MODEL>      Observations_;
  typedef oops::ObsOperator<MODEL>       ObsOperator_;
  typedef oops::ObsSpaces<MODEL>         ObsSpaces_;
  typedef oops::ObsVector<MODEL>         ObsVector_;

 public:
// -----------------------------------------------------------------------------
  explicit RunCRTM(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {}
// -----------------------------------------------------------------------------
  virtual ~RunCRTM() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {
//  Setup observation window
    const util::DateTime winbgn(fullConfig.getString("window begin"));
    const util::Duration winlen(fullConfig.getString("window length"));
    const util::DateTime winend(winbgn + winlen);

//  Setup observations
    ObsSpaces_ obsdb(fullConfig, this->getComm(), winbgn, winend);

    oops::Variables diagvars;

    std::vector<eckit::LocalConfiguration> conf;
    fullConfig.get("observations", conf);

    for (std::size_t jj = 0; jj < obsdb.size(); ++jj) {
      eckit::LocalConfiguration obsopconf(conf[jj], "obs operator");

      ObsOperator_ hop(obsdb[jj], obsopconf);

      const eckit::LocalConfiguration gconf(conf[jj], "geovals");
      const GeoVaLs_ gval(gconf, obsdb[jj], hop.requiredVars());

      const ObsAuxCtrl_ ybias(obsdb[jj], conf[jj]);

      ObsVector_ hofx(obsdb[jj]);
      ObsDiags_ diag(obsdb[jj], hop.locations(), diagvars);

      hop.simulateObs(gval, hofx, ybias, diag);

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
    return "oops::RunCRTM<" + MODEL::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace ufo

#endif  // MAINS_RUNCRTM_H_
