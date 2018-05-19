/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

//TODO: through the file replace UFO_OBSEXAMPLETLAD_H with the unique string (e.g. UFO_OBS<YOUR_OBS_OPERATOR_NAME>TLAD_H
#ifndef UFO_OBSEXAMPLETLAD_H_
#define UFO_OBSEXAMPLETLAD_H_

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
/// Example for observation operator TL and AD class
//  TODO: through the file replace ObsExampleTLAD with Obs<Your_Obs_Operator_Name>TLAD
template <typename MODEL>
class ObsExampleTLAD : public oops::LinearObsOperBase<MODEL>,
                        private util::ObjectCounter<ObsExampleTLAD<MODEL>> {
public:
  static const std::string classname() {return "ufo::ObsExampleTLAD";}

  ObsExampleTLAD(const ioda::ObsSpace &, const eckit::Configuration &);
  virtual ~ObsExampleTLAD();

  // Obs Operators
  void setTrajectory(const GeoVaLs &, const ObsBias &);
  void obsEquivTL(const GeoVaLs &, ioda::ObsVector &, const ObsBiasIncrement &) const;
  void obsEquivAD(GeoVaLs &, const ioda::ObsVector &, ObsBiasIncrement &) const;

  // Other
  const oops::Variables & variables() const {return *varin_;}

  int & toFortran() {return keyOper_;}
  const int & toFortran() const {return keyOper_;}

private:
  void print(std::ostream &) const;
  F90hop keyOper_;
  const ioda::ObsSpace& odb_;
  boost::scoped_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------
template <typename MODEL>
ObsExampleTLAD<MODEL>::ObsExampleTLAD(const ioda::ObsSpace & odb, const eckit::Configuration & config)
  : keyOper_(0), varin_(), odb_(odb)
{
  const eckit::Configuration * configc = &config;
  // TODO: replace ufo_example_tlad_setup_f90 with the call to your Fortran routine
  //       to setup tl/ad obs operator (defined in ObsExampleTLAD.interface.F90)
  ufo_example_tlad_setup_f90(keyOper_, &configc);
  // TODO: list the variables for GeoVaLs that are needed for the observation 
  //       operator below in vv (e.g., vv{"temperature", "humidity"})
  const std::vector<std::string> vv{""};
  varin_.reset(new oops::Variables(vv));
  oops::Log::trace() << "ObsExampleTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL>
ObsExampleTLAD<MODEL>::~ObsExampleTLAD() {
  // TODO: replace ufo_example_tlad_delete_f90 with the call to your Fortran routine
  //       to destruct tl/ad observation operator (defined in ObsExampleTLAD.interface.F90)
  ufo_example_tlad_delete_f90(keyOper_);
  oops::Log::trace() << "ObsExampleTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsExampleTLAD<MODEL>::setTrajectory(const GeoVaLs & geovals, const ObsBias & bias) {
  // TODO: replace ufo_example_tlad_settraj_f90 with the call to your Fortran routine
  //       to set trajectory for tl/ad (defined in ObsExampleTLAD.interface.F90)
  ufo_example_tlad_settraj_f90(keyOper_, geovals.toFortran(), odb_.toFortran());
  oops::Log::trace() << "ObsExampleTLAD: trajectory set" << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsExampleTLAD<MODEL>::obsEquivTL(const GeoVaLs & geovals, ioda::ObsVector & ovec,
                             const ObsBiasIncrement & bias) const {
  // TODO: replace ufo_example_tlad_eqv_tl_f90 with the call to your Fortran routine
  //       to apply tl observation operator (defined in ObsExampleTLAD.interface.F90)
  ufo_example_tlad_eqv_tl_f90(keyOper_, geovals.toFortran(), odb_.toFortran(), ovec.toFortran());
  oops::Log::trace() << "ObsExampleTLAD: tangent linear observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsExampleTLAD<MODEL>::obsEquivAD(GeoVaLs & geovals, const ioda::ObsVector & ovec,
                             ObsBiasIncrement & bias) const {
  // TODO: replace ufo_example_tlad_eqv_ad_f90 with the call to your Fortran routine
  //       to apply ad observation operator (defined in ObsExampleTLAD.interface.F90)
  ufo_example_tlad_eqv_ad_f90(keyOper_, geovals.toFortran(), odb_.toFortran(), ovec.toFortran());
  oops::Log::trace() << "ObsExampleTLAD: adjoint observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsExampleTLAD<MODEL>::print(std::ostream & os) const {
  os << "ObsExampleTLAD::print not implemented" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OBSEXAMPLETLAD_H_
