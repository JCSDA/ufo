/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ObsStericHeightTLAD.h"

#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "oops/base/Variables.h"
#include "ufo/LinearObsOperatorBase.h"
#include "ioda/ObsSpace.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Logger.h"
#include "ioda/ObsVector.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsBiasIncrement.h"

namespace ufo {

  // -----------------------------------------------------------------------------
static LinearObsOperatorMaker<ObsStericHeightTLAD> makerObsStericHeightTLAD_("ObsStericHeightTLAD");
  // -----------------------------------------------------------------------------

    ObsStericHeightTLAD::ObsStericHeightTLAD(const ioda::ObsSpace & odb, const eckit::Configuration & config)
    : keyOperStericHeight_(0), varin_(), traj_()
    {
      std::cout << "steric height tlad =============================" << std::endl;
      const eckit::Configuration * configc = &config;
      ufo_stericheight_tlad_setup_f90(keyOperStericHeight_, &configc);
      const std::vector<std::string> vv{"sea_surface_height_above_geoid",
	  "ocean_potential_temperature",
	  "ocean_salinity"};
      varin_.reset(new oops::Variables(vv));
      traj_.reset(new GeoVaLs(config, oops::Variables(vv)));

      oops::Variables vars(vv);
      GeoVaLs traj(config,vars);  
      //ufo_stericheight_tlad_gettraj_f90(keyOperStericHeight_, odb.nobs(), vars.toFortran(), traj.toFortran());
      oops::Log::trace() << "ObsStericHeightTLAD created" << std::endl;
    }

  // -----------------------------------------------------------------------------

    ObsStericHeightTLAD::~ObsStericHeightTLAD() {
    ufo_stericheight_tlad_delete_f90(keyOperStericHeight_);
    oops::Log::trace() << "ObsStericHeightTLAD destrcuted" << std::endl;  
  }

  // -----------------------------------------------------------------------------

    void ObsStericHeightTLAD::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias) {
    std::cout << "steric height tlad settraj =============================" << std::endl;  
    ufo_stericheight_tlad_settraj_f90(keyOperStericHeight_, geovals.toFortran());
    oops::Log::trace() << "ObsStericHeightTLAD trajectory was set " << geovals << std::endl;  
  }

  // -----------------------------------------------------------------------------

    void ObsStericHeightTLAD::simulateObsTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
						const ObsBiasIncrement & bias) const {
    ufo_stericheight_tlad_eqv_tl_f90(keyOperStericHeight_, geovals.toFortran(), ovec.toFortran());
  }

  // -----------------------------------------------------------------------------

    void ObsStericHeightTLAD::simulateObsAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
						ObsBiasIncrement & bias) const {
    ufo_stericheight_tlad_eqv_ad_f90(keyOperStericHeight_, geovals.toFortran(), ovec.toFortran());
  }

  // -----------------------------------------------------------------------------

    void ObsStericHeightTLAD::print(std::ostream & os) const {
    os << "ObsStericHeightTLAD::print not implemented" << std::endl;
  }

  // -----------------------------------------------------------------------------

}  // namespace ufo
