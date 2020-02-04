/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_INSTANTIATEOBSFILTERFACTORY_H_
#define UFO_INSTANTIATEOBSFILTERFACTORY_H_

#include "oops/base/instantiateObsFilterFactory.h"
#include "oops/interface/ObsFilter.h"
#include "ufo/filters/BackgroundCheck.h"
#include "ufo/filters/BlackList.h"
#include "ufo/filters/DifferenceCheck.h"
#include "ufo/filters/Gaussian_Thinning.h"
#include "ufo/filters/MWCLWCheck.h"
#include "ufo/filters/ObsBoundsCheck.h"
#include "ufo/filters/ObsDomainCheck.h"
#include "ufo/filters/ObsDomainErrCheck.h"
#include "ufo/filters/PreQC.h"
#include "ufo/filters/QCmanager.h"
#include "ufo/filters/Thinning.h"
#include "ufo/gnssro/QC/BackgroundCheckRONBAM.h"
#include "ufo/gnssro/QC/ROobserror.h"
#include "ufo/onedvarcheck/OneDVarCheck.h"

namespace ufo {
template<typename MODEL> void instantiateObsFilterFactory() {
  oops::instantiateObsFilterFactory<MODEL>();
  static oops::FilterMaker<MODEL, oops::ObsFilter<MODEL, ufo::QCmanager> >
           makerChk0_("QCmanager");
  static oops::FilterMaker<MODEL, oops::ObsFilter<MODEL, ufo::PreQC> >
           makerChk1_("PreQC");
  static oops::FilterMaker<MODEL, oops::ObsFilter<MODEL, ufo::ObsDomainCheck> >
           makerChk2_("Domain Check");
  static oops::FilterMaker<MODEL, oops::ObsFilter<MODEL, ufo::ObsBoundsCheck> >
           makerChk3_("Bounds Check");
  static oops::FilterMaker<MODEL, oops::ObsFilter<MODEL, ufo::BlackList> >
           makerChk4_("BlackList");
  static oops::FilterMaker<MODEL, oops::ObsFilter<MODEL, ufo::BackgroundCheck> >
           makerChk5_("Background Check");
  static oops::FilterMaker<MODEL, oops::ObsFilter<MODEL, ufo::DifferenceCheck> >
           makerChk6_("Difference Check");
  static oops::FilterMaker<MODEL, oops::ObsFilter<MODEL, ufo::ROobserror> >
           makerChk7_("ROobserror");
  static oops::FilterMaker<MODEL, oops::ObsFilter<MODEL, ufo::Thinning> >
           makerChk8_("Thinning");
  static oops::FilterMaker<MODEL, oops::ObsFilter<MODEL, ufo::Gaussian_Thinning> >
           makerChk9_("Gaussian_Thinning");
  static oops::FilterMaker<MODEL, oops::ObsFilter<MODEL, ufo::MWCLWCheck> >
           makerChk10_("MWCLW Check");
  static oops::FilterMaker<MODEL, oops::ObsFilter<MODEL, ufo::ObsDomainErrCheck> >
           makerChk11_("DomainErr Check");
  static oops::FilterMaker<MODEL, oops::ObsFilter<MODEL, ufo::BackgroundCheckRONBAM> >
           makerChk13_("Background Check RONBAM");
  static oops::FilterMaker<MODEL, oops::ObsFilter<MODEL, ufo::OneDVarCheck> >
           makerChk14_("OneDVar Check");
}

}  // namespace ufo

#endif  // UFO_INSTANTIATEOBSFILTERFACTORY_H_
