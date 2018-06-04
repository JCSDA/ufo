/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OBSADT_H_
#define UFO_OBSADT_H_

#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "eckit/config/Configuration.h"
#include "oops/base/Variables.h"
#include "oops/interface/ObsOperatorBase.h"
#include "ioda/ObsSpace.h"
#include "ufo/GeoVaLs.h"
#include "ioda/Locations.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsBiasIncrement.h"
#include "ioda/ObsVector.h"
#include "oops/util/ObjectCounter.h"

namespace ufo {

// -----------------------------------------------------------------------------
/// adt observation for UFO.
template <typename MODEL>
class ObsADT : public oops::ObsOperatorBase<MODEL>,
                  private util::ObjectCounter<ObsADT<MODEL>> {
 public:
  static const std::string classname() {return "ufo::ObsADT";}

  ObsADT(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsADT();

// Obs Operator
  void obsEquiv(const GeoVaLs &, ioda::ObsVector &, const ObsBias &) const;

// Other
  const oops::Variables & variables() const {return *varin_;}

  int & toFortran() {return keyOperADT_;}
  const int & toFortran() const {return keyOperADT_;}

 private:
  void print(std::ostream &) const;
  F90hop keyOperADT_;
  const ioda::ObsSpace& odb_;
  boost::scoped_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------
template <typename MODEL>
  ObsADT<MODEL>::ObsADT(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : keyOperADT_(0), varin_(), odb_(odb)
{
  const eckit::Configuration * configc = &config;
  ufo_adt_setup_f90(keyOperADT_, &configc);
  const std::vector<std::string> vv{"sea_surface_height_above_geoid"};
  varin_.reset(new oops::Variables(vv));
  oops::Log::trace() << "ObsADT created." << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL>
ObsADT<MODEL>::~ObsADT() {
  ufo_adt_delete_f90(keyOperADT_);
  oops::Log::trace() << "ObsADT destructed" << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL>
  void ObsADT<MODEL>::obsEquiv(const GeoVaLs & gom, ioda::ObsVector & ovec,
                             const ObsBias & bias) const {
  ufo_adt_eqv_f90(keyOperADT_, gom.toFortran(), odb_.toFortran(), ovec.toFortran(), bias.toFortran());
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsADT<MODEL>::print(std::ostream & os) const {
  os << "ObsADT::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OBSADT_H_
