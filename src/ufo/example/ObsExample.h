/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */


//TODO: through the file replace UFO_OBSEXAMPLE_H with the unique string (e.g. UFO_OBS<YOUR_OBS_OPERATOR_NAME>_H
#ifndef UFO_OBSEXAMPLE_H_
#define UFO_OBSEXAMPLE_H_

#include <ostream>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "eckit/config/Configuration.h"
#include "oops/base/Variables.h"
#include "oops/interface/ObsOperatorBase.h"
#include "ufo/ObsSpace.h"
#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsBiasIncrement.h"
#include "ufo/ObsVector.h"
#include "util/ObjectCounter.h"

namespace ufo {

// -----------------------------------------------------------------------------
/// Example for the observation operator class.
//  TODO: through the file replace ObsExample with Obs<Your_Obs_Operator_Name>
template <typename MODEL>
class ObsExample : public oops::ObsOperatorBase<MODEL>,
                  private util::ObjectCounter<ObsExample<MODEL>> {
 public:
  static const std::string classname() {return "ufo::ObsExample";}

  ObsExample(const ObsSpace &, const eckit::Configuration &);
  virtual ~ObsExample();

  // Obs Operator
  void obsEquiv(const GeoVaLs &, ObsVector &, const ObsBias &) const;

  // Other
  const oops::Variables & variables() const {return *varin_;}

  int & toFortran() {return keyOper_;}
  const int & toFortran() const {return keyOper_;}

 private:
  void print(std::ostream &) const;
  F90hop keyOper_;
  const ObsSpace& odb_;
  boost::scoped_ptr<const oops::Variables> varin_;
};

// -----------------------------------------------------------------------------
template <typename MODEL>
ObsExample<MODEL>::ObsExample(const ObsSpace & odb, const eckit::Configuration & config)
  : keyOper_(0), varin_(), odb_(odb)
{
  const eckit::Configuration * configc = &config;
  // TODO: replace ufo_example_setup_f90 with the call to your Fortran routine
  //       to setup obs operator (defined in ObsExample.interface.F90)
  ufo_example_setup_f90(keyOper_, &configc);
  // TODO: list the variables for GeoVaLs that are needed for the observation 
  //       operator below in vv (e.g., vv{"temperature", "humidity"})
  const std::vector<std::string> vv{""};
  varin_.reset(new oops::Variables(vv));
  oops::Log::trace() << "ObsExample created." << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL>
ObsExample<MODEL>::~ObsExample() {
  // TODO: replace ufo_example_delete_f90 with the call to your Fortran routine
  //       to destruct observation operator (defined in ObsExample.interface.F90)
  ufo_example_delete_f90(keyOper_);
  oops::Log::trace() << "ObsExample destructed" << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsExample<MODEL>::obsEquiv(const GeoVaLs & gv, ObsVector & ovec,
                                  const ObsBias & bias) const {
  // TODO: replace ufo_example_eqv_f90 with the call to your Fortran routine
  //       to apply observation operator (defined in ObsExample.interface.F90)

  ufo_example_eqv_f90(keyOper_, gv.toFortran(), odb_.toFortran(), ovec.toFortran(), bias.toFortran());
  oops::Log::trace() << "ObsExample: observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL>
void ObsExample<MODEL>::print(std::ostream & os) const {
  os << "ObsExample::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
#endif  // UFO_OBSEXAMPLE_H_
