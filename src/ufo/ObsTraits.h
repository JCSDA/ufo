/*
 * (C) Copyright 2017-2020 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OBSTRAITS_H_
#define UFO_OBSTRAITS_H_

#include <string>

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "AnalyticInit.h"
#include "Locations.h"
#include "GeoVaLs.h"
#include "ObsBias.h"
#include "ObsBiasCovariance.h"
#include "ObsBiasIncrement.h"
#include "ObsDiagnostics.h"
#include "ObsOperator.h"
#include "LinearObsOperator.h"

namespace ufo {

struct ObsTraits {
  static std::string name() {return "UFO and IODA observations";}

  typedef ufo::AnalyticInit        AnalyticInit;
  typedef ufo::GeoVaLs             GeoVaLs;
  typedef ufo::ObsDiagnostics      ObsDiagnostics;
  typedef ufo::Locations           Locations;
  typedef ioda::ObsSpace           ObsSpace;
  typedef ioda::ObsVector          ObsVector;
  template <typename DATATYPE> using ObsDataVector = ioda::ObsDataVector<DATATYPE>;

  typedef ufo::ObsOperator         ObsOperator;
  typedef ufo::LinearObsOperator   LinearObsOperator;

  typedef ufo::ObsBias             ObsAuxControl;
  typedef ufo::ObsBiasIncrement    ObsAuxIncrement;
  typedef ufo::ObsBiasCovariance   ObsAuxCovariance;
};

}  // namespace ufo

#endif  // UFO_OBSTRAITS_H_
