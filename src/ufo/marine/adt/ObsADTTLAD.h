/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_OBSADTTLAD_H_
#define UFO_OBSADTTLAD_H_

#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "oops/base/Variables.h"
#include "oops/interface/LinearObsOperBase.h"
#include "ioda/ObsSpace.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Logger.h"

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
/// ADT observation for  model.
template <typename MODEL>
class ObsADTTLAD : public oops::LinearObsOperBase<MODEL>, 
                              private util::ObjectCounter<ObsADTTLAD<MODEL>> {
public:
  static const std::string classname() {return "ufo::ObsADTTLAD";}

  ObsADTTLAD(const ioda::ObsSpace &, const eckit::Configuration &);    
  virtual ~ObsADTTLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, const ObsBias &);
  void obsEquivTL(const GeoVaLs &, ioda::ObsVector &, const ObsBiasIncrement &) const;
  void obsEquivAD(GeoVaLs &, const ioda::ObsVector &, ObsBiasIncrement &) const;

  // Other
  const oops::Variables & variables() const {return *varin_;}

  int & toFortran() {return keyOperADT_;}
  const int & toFortran() const {return keyOperADT_;}

private:
  void print(std::ostream &) const;
  F90hop keyOperADT_;
  boost::scoped_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------
template <typename MODEL>
  ObsADTTLAD<MODEL>::ObsADTTLAD(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : keyOperADT_(0), varin_()
{
  const eckit::Configuration * configc = &config;
  ufo_adt_tlad_setup_f90(keyOperADT_, &configc);
  const std::vector<std::string> vv{"sea_surface_height_above_geoid"};
  varin_.reset(new oops::Variables(vv));
  oops::Log::trace() << "ObsADTTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL>
ObsADTTLAD<MODEL>::~ObsADTTLAD() {
  ufo_adt_tlad_delete_f90(keyOperADT_);
  oops::Log::trace() << "ObsADTTLAD destrcuted" << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsADTTLAD<MODEL>::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias) {
  ufo_adt_tlad_settraj_f90(keyOperADT_, geovals.toFortran());
}

// -----------------------------------------------------------------------------
template <typename MODEL>
  void ObsADTTLAD<MODEL>::obsEquivTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                               const ObsBiasIncrement & bias) const {
  ufo_adt_tlad_eqv_tl_f90(keyOperADT_, geovals.toFortran(), ovec.toFortran());
}

// -----------------------------------------------------------------------------
template <typename MODEL>
  void ObsADTTLAD<MODEL>::obsEquivAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                               ObsBiasIncrement & bias) const {
  ufo_adt_tlad_eqv_ad_f90(keyOperADT_, geovals.toFortran(), ovec.toFortran());
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsADTTLAD<MODEL>::print(std::ostream & os) const {
  os << "ObsADTTLAD::print not implemented" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OBSADTTLAD_H_
