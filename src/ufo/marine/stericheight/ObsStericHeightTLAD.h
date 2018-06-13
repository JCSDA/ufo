/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSSTERICHEIGHTTLAD_H_
#define UFO_OBSSTERICHEIGHTTLAD_H_

#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "oops/base/Variables.h"
#include "oops/interface/LinearObsOperBase.h"
#include "ioda/ObsSpace.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Logger.h"
#include "ufo/FortranMarine.h"

// Forward declarations
namespace util {
  class DateTime;
}

namespace ioda {
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsBias;
  class ObsBiasIncrement;

  // -----------------------------------------------------------------------------
  /// Simulated Steric height for  model.
  template <typename MODEL>
    class ObsStericHeightTLAD : public oops::LinearObsOperBase<MODEL>, 
    private util::ObjectCounter<ObsStericHeightTLAD<MODEL>> {
  public:
      static const std::string classname() {return "ufo::ObsStericHeightTLAD";}

      ObsStericHeightTLAD(const ioda::ObsSpace &, const eckit::Configuration &);    
      virtual ~ObsStericHeightTLAD();

      // Obs Operators
      void setTrajectory(const GeoVaLs &, const ObsBias &);
      void obsEquivTL(const GeoVaLs &, ioda::ObsVector &, const ObsBiasIncrement &) const;
      void obsEquivAD(GeoVaLs &, const ioda::ObsVector &, ObsBiasIncrement &) const;

      // Other
      const oops::Variables & variables() const {return *varin_;}
      int & toFortran() {return keyOperStericHeight_;}
      const int & toFortran() const {return keyOperStericHeight_;}
      //const GeoVaLs * traj_;
  
  private:
      void print(std::ostream &) const;
      F90hop keyOperStericHeight_;
      boost::scoped_ptr<const GeoVaLs> traj_;
      boost::scoped_ptr<const oops::Variables> varin_;
    };

  // -----------------------------------------------------------------------------
  template <typename MODEL>
    ObsStericHeightTLAD<MODEL>::ObsStericHeightTLAD(const ioda::ObsSpace & odb, const eckit::Configuration & config)
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
  template <typename MODEL>
    ObsStericHeightTLAD<MODEL>::~ObsStericHeightTLAD() {
    ufo_stericheight_tlad_delete_f90(keyOperStericHeight_);
    oops::Log::trace() << "ObsStericHeightTLAD destrcuted" << std::endl;  
  }

  // -----------------------------------------------------------------------------
  template <typename MODEL>
    void ObsStericHeightTLAD<MODEL>::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias) {
    std::cout << "steric height tlad settraj =============================" << std::endl;  
    ufo_stericheight_tlad_settraj_f90(keyOperStericHeight_, geovals.toFortran());
    oops::Log::trace() << "ObsStericHeightTLAD trajectory was set " << geovals << std::endl;  
  }

  // -----------------------------------------------------------------------------
  template <typename MODEL>
    void ObsStericHeightTLAD<MODEL>::obsEquivTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
						const ObsBiasIncrement & bias) const {
    ufo_stericheight_tlad_eqv_tl_f90(keyOperStericHeight_, geovals.toFortran(), ovec.toFortran());
  }

  // -----------------------------------------------------------------------------
  template <typename MODEL>
    void ObsStericHeightTLAD<MODEL>::obsEquivAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
						ObsBiasIncrement & bias) const {
    ufo_stericheight_tlad_eqv_ad_f90(keyOperStericHeight_, geovals.toFortran(), ovec.toFortran());
  }

  // -----------------------------------------------------------------------------
  template <typename MODEL>
    void ObsStericHeightTLAD<MODEL>::print(std::ostream & os) const {
    os << "ObsStericHeightTLAD::print not implemented" << std::endl;
  }
  // -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OBSSTERICHEIGHTTLAD_H_
