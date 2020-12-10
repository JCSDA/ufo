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
#include "ufo/filters/MetOfficeBuddyCheck.h"
#include "ufo/filters/MWCLWCheck.h"
#include "ufo/filters/ObsBoundsCheck.h"
#include "ufo/filters/ObsDerivativeCheck.h"
#include "ufo/filters/ObsDiagnosticsWriter.h"
#include "ufo/filters/ObsDomainCheck.h"
#include "ufo/filters/ObsDomainErrCheck.h"
#include "ufo/filters/PoissonDiskThinning.h"
#include "ufo/filters/PreQC.h"
#include "ufo/filters/ProfileConsistencyChecks.h"
#include "ufo/filters/QCmanager.h"
#include "ufo/filters/RatioCheck.h"
#include "ufo/filters/TemporalThinning.h"
#include "ufo/filters/Thinning.h"
#include "ufo/filters/TrackCheck.h"
#include "ufo/filters/TrackCheckShip.h"
#include "ufo/filters/variabletransforms/WindComponents.h"
#include "ufo/filters/variabletransforms/WindSpeedAndDirection.h"
#include "ufo/gnssro/QC/BackgroundCheckRONBAM.h"
#include "ufo/gnssro/QC/ROobserror.h"

#if defined(RTTOV_FOUND)
  #include "ufo/filters/rttovonedvarcheck/RTTOVOneDVarCheck.h"
#endif

namespace ufo {
template<typename MODEL> void instantiateObsFilterFactory() {
  oops::instantiateObsFilterFactory<MODEL>();
  static oops::FilterMaker<MODEL, oops::ObsFilter<MODEL, ufo::QCmanager> >
           qcManagerMaker("QCmanager");
  static oops::FilterMaker<MODEL, oops::ObsFilter<MODEL, ufo::PreQC> >
           preQCMaker("PreQC");
  static oops::FilterMaker<MODEL, oops::ObsFilter<MODEL, ufo::ObsDomainCheck> >
           domainCheckMaker("Domain Check");
  static oops::FilterMaker<MODEL, oops::ObsFilter<MODEL, ufo::ObsBoundsCheck> >
           boundsCheckMaker("Bounds Check");
  static oops::FilterMaker<MODEL, oops::ObsFilter<MODEL, ufo::BlackList> >
           blackListMaker("BlackList");
  static oops::FilterMaker<MODEL, oops::ObsFilter<MODEL, ufo::BackgroundCheck> >
           backgroundCheckMaker("Background Check");
  static oops::FilterMaker<MODEL, oops::ObsFilter<MODEL, ufo::RatioCheck> >
           ratioCheckMaker("Ratio Check");
  static oops::FilterMaker<MODEL, oops::ObsFilter<MODEL, ufo::DifferenceCheck> >
           differenceCheckMaker("Difference Check");
  static oops::FilterMaker<MODEL, oops::ObsFilter<MODEL, ufo::ROobserror> >
           ROobserrorMaker("ROobserror");
  static oops::FilterMaker<MODEL, oops::ObsFilter<MODEL, ufo::Thinning> >
           thinningMaker("Thinning");
  static oops::FilterMaker<MODEL, oops::ObsFilter<MODEL, ufo::Gaussian_Thinning> >
           gaussianThinningMaker("Gaussian Thinning");
  static oops::FilterMaker<MODEL, oops::ObsFilter<MODEL, ufo::MWCLWCheck> >
           MWCLWCheckMaker("MWCLW Check");
  static oops::FilterMaker<MODEL, oops::ObsFilter<MODEL, ufo::ObsDomainErrCheck> >
           domainErrCheckMaker("DomainErr Check");
  static oops::FilterMaker<MODEL, oops::ObsFilter<MODEL, ufo::ProfileConsistencyChecks> >
           profileConsistencyChecksMaker("Profile Consistency Checks");
  static oops::FilterMaker<MODEL, oops::ObsFilter<MODEL, ufo::BackgroundCheckRONBAM> >
           backgroundCheckRONBAMMaker("Background Check RONBAM");
  static oops::FilterMaker<MODEL, oops::ObsFilter<MODEL, ufo::TemporalThinning> >
           temporalThinningMaker("Temporal Thinning");
  static oops::FilterMaker<MODEL, oops::ObsFilter<MODEL, ufo::PoissonDiskThinning> >
           poissonDiskThinningMaker("Poisson Disk Thinning");
  static oops::FilterMaker<MODEL, oops::ObsFilter<MODEL, ufo::ObsDiagnosticsWriter> >
           YDIAGsaverMaker("YDIAGsaver");
  static oops::FilterMaker<MODEL, oops::ObsFilter<MODEL, ufo::TrackCheck> >
           TrackCheckMaker("Track Check");
  static oops::FilterMaker<MODEL, oops::ObsFilter<MODEL, ufo::MetOfficeBuddyCheck> >
           MetOfficeBuddyCheckMaker("Met Office Buddy Check");
  static oops::FilterMaker<MODEL, oops::ObsFilter<MODEL, ufo::ObsDerivativeCheck> >
           DerivativeCheckMaker("Derivative Check");
  static oops::FilterMaker<MODEL, oops::ObsFilter<MODEL, ufo::TrackCheckShip> >
           ShipTrackCheckMaker("Ship Track Check");
  static oops::FilterMaker<MODEL, oops::ObsFilter<MODEL, ufo::WindComponents> >
           WindComponentsMaker("Wind Components");
  static oops::FilterMaker<MODEL, oops::ObsFilter<MODEL, ufo::WindSpeedAndDirection> >
           WindSpeedAndDirectionMaker("Wind Speed And Direction");

  // Only include this filter if rttov is present
  #if defined(RTTOV_FOUND)
    static oops::FilterMaker<MODEL, oops::ObsFilter<MODEL, ufo::RTTOVOneDVarCheck> >
             RTTOVOneDVarCheckMaker("RTTOV OneDVar Check");
  #endif

  // For backward compatibility, register some filters under legacy names used in the past
  static oops::FilterMaker<MODEL, oops::ObsFilter<MODEL, ufo::Gaussian_Thinning> >
           legacyGaussianThinningMaker("Gaussian_Thinning");
  static oops::FilterMaker<MODEL, oops::ObsFilter<MODEL, ufo::TemporalThinning> >
           legacyTemporalThinningMaker("TemporalThinning");
}

}  // namespace ufo

#endif  // UFO_INSTANTIATEOBSFILTERFACTORY_H_
