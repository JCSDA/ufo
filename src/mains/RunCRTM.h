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
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace ufo {

template <typename MODEL> class RunCRTM : public oops::Application {
  typedef oops::GeoVaLs<MODEL>           GeoVaLs_;
  typedef oops::ObsAuxControl<MODEL>     ObsAuxCtrl_;
  typedef oops::Observations<MODEL>      Observations_;
  typedef oops::ObsOperator<MODEL>       ObsOperator_;
  typedef oops::ObsSpaces<MODEL>         ObsSpaces_;
  typedef oops::ObsVector<MODEL>         ObsVector_;

 public:
// -----------------------------------------------------------------------------
  RunCRTM() {}
// -----------------------------------------------------------------------------
  virtual ~RunCRTM() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {
//  Setup observation window
    const eckit::LocalConfiguration windowConf(fullConfig, "Assimilation Window");
    const util::Duration winlen(windowConf.getString("Length"));
    const util::DateTime winbgn(windowConf.getString("Begin"));
    const util::DateTime winend(winbgn + winlen);
    oops::Log::info() << "Observation window is:" << windowConf << std::endl;

//  Setup observations
    eckit::LocalConfiguration obsconf(fullConfig, "Observations");
    oops::Log::debug() << "Observations configuration is:" << obsconf << std::endl;
    ObsSpaces_ obsdb(obsconf, winbgn, winend);

    std::vector<eckit::LocalConfiguration> conf;
    obsconf.get("ObsTypes", conf);

    for (std::size_t jj = 0; jj < obsdb.size(); ++jj) {
      ObsOperator_ hop(obsdb[jj]);

      const eckit::LocalConfiguration gconf(conf[jj], "GeoVaLs");
      const GeoVaLs_ gval(gconf, hop.variables());

      eckit::LocalConfiguration biasConf;
      conf[jj].get("ObsBias", biasConf);
      const ObsAuxCtrl_ ybias(biasConf);

      ObsVector_ ovec(obsdb[jj], hop.observed());

      hop.simulateObs(gval, ovec, ybias);

      const double zz = ovec.rms();
      const double xx = conf[jj].getDouble("rmsequiv");
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
