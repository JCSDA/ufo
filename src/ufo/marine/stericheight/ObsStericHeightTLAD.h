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
#include "ufo/ObsSpace.h"
#include "util/ObjectCounter.h"
#include "util/Logger.h"

// Forward declarations
namespace util {
  class DateTime;
}

namespace ufo {
  class GeoVaLs;
  class ObsBias;
  class ObsBiasIncrement;
  class ObsVector;

  // -----------------------------------------------------------------------------
  /// Simulated Steric height for  model.
  template <typename MODEL>
    class ObsStericHeightTLAD : public oops::LinearObsOperBase<MODEL>, 
    private util::ObjectCounter<ObsStericHeightTLAD<MODEL>> {
  public:
      static const std::string classname() {return "ufo::ObsStericHeightTLAD";}

      ObsStericHeightTLAD(const ObsSpace &, const eckit::Configuration &);    
      virtual ~ObsStericHeightTLAD();

      // Obs Operators
      void setTrajectory(const GeoVaLs &, const ObsBias &);
      void obsEquivTL(const GeoVaLs &, ObsVector &, const ObsBiasIncrement &) const;
      void obsEquivAD(GeoVaLs &, const ObsVector &, ObsBiasIncrement &) const;

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
    ObsStericHeightTLAD<MODEL>::ObsStericHeightTLAD(const ObsSpace & odb, const eckit::Configuration & config)
    : keyOperStericHeight_(0), varin_()
    {
      std::cout << "steric height tlad =============================" << std::endl;
      const eckit::Configuration * configc = &config;
      ufo_stericheight_tlad_setup_f90(keyOperStericHeight_, &configc);
      const std::vector<std::string> vv{"sea_surface_height_above_geoid",
	  "ocean_potential_temperature",
	  "ocean_salinity"};
      varin_.reset(new oops::Variables(vv));
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
    ufo_stericheight_tlad_settraj_f90(keyOperStericHeight_, geovals.toFortran());
    oops::Log::trace() << "ObsStericHeightTLAD trajectory was set " << geovals << std::endl;  
  }

  // -----------------------------------------------------------------------------
  template <typename MODEL>
    void ObsStericHeightTLAD<MODEL>::obsEquivTL(const GeoVaLs & geovals, ObsVector & ovec,
						const ObsBiasIncrement & bias) const {
    ufo_stericheight_tlad_eqv_tl_f90(keyOperStericHeight_, geovals.toFortran(), ovec.toFortran());
  }

  // -----------------------------------------------------------------------------
  template <typename MODEL>
    void ObsStericHeightTLAD<MODEL>::obsEquivAD(GeoVaLs & geovals, const ObsVector & ovec,
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
