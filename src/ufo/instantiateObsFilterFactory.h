/*
 * (C) Copyright 2019-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_INSTANTIATEOBSFILTERFACTORY_H_
#define UFO_INSTANTIATEOBSFILTERFACTORY_H_

#include "oops/base/instantiateObsFilterFactory.h"
#include "oops/interface/ObsFilter.h"
#include "ufo/filters/AcceptList.h"
#include "ufo/filters/BackgroundCheck.h"
#include "ufo/filters/BayesianBackgroundCheck.h"
#include "ufo/filters/BayesianBackgroundQCFlags.h"
#include "ufo/filters/BlackList.h"
#include "ufo/filters/DifferenceCheck.h"
#include "ufo/filters/FinalCheck.h"
#include "ufo/filters/Gaussian_Thinning.h"
#include "ufo/filters/gnssroonedvarcheck/GNSSROOneDVarCheck.h"
#include "ufo/filters/HistoryCheck.h"
#include "ufo/filters/ImpactHeightCheck.h"
#include "ufo/filters/MetOfficeBuddyCheck.h"
#include "ufo/filters/ModelObThreshold.h"
#include "ufo/filters/MWCLWCheck.h"
#include "ufo/filters/ObsBoundsCheck.h"
#include "ufo/filters/ObsDerivativeCheck.h"
#include "ufo/filters/ObsDiagnosticsWriter.h"
#include "ufo/filters/ObsDomainCheck.h"
#include "ufo/filters/ObsDomainErrCheck.h"
#include "ufo/filters/PerformAction.h"
#include "ufo/filters/PoissonDiskThinning.h"
#include "ufo/filters/PreQC.h"
#include "ufo/filters/ProbabilityGrossErrorWholeReport.h"
#include "ufo/filters/ProfileBackgroundCheck.h"
#include "ufo/filters/ProfileConsistencyChecks.h"
#include "ufo/filters/ProfileFewObsCheck.h"
#include "ufo/filters/QCmanager.h"
#include "ufo/filters/SatName.h"
#include "ufo/filters/SatwindInversionCorrection.h"
#include "ufo/filters/StuckCheck.h"
#include "ufo/filters/TemporalThinning.h"
#include "ufo/filters/Thinning.h"
#include "ufo/filters/TrackCheck.h"
#include "ufo/filters/TrackCheckShip.h"
#include "ufo/filters/VariableAssignment.h"
#include "ufo/filters/VariableTransforms.h"
#include "ufo/gnssro/QC/BackgroundCheckRONBAM.h"
#include "ufo/gnssro/QC/ROobserror.h"

#if defined(RTTOV_FOUND)
  #include "ufo/filters/rttovonedvarcheck/RTTOVOneDVarCheck.h"
#endif

namespace ufo {
template<typename OBS> void instantiateObsFilterFactory() {
  oops::instantiateObsFilterFactory<OBS>();
  static oops::FilterMaker<OBS, oops::ObsFilter<OBS, ufo::QCmanager> >
           qcManagerMaker("QCmanager");
  static oops::FilterMaker<OBS, oops::ObsFilter<OBS, ufo::FinalCheck> >
           finalCheckMaker("Final Check");
  static oops::FilterMaker<OBS, oops::ObsFilter<OBS, ufo::PreQC> >
           preQCMaker("PreQC");
  static oops::FilterMaker<OBS, oops::ObsFilter<OBS, ufo::ObsDomainCheck> >
           domainCheckMaker("Domain Check");
  static oops::FilterMaker<OBS, oops::ObsFilter<OBS, ufo::SatName> >
           satnameCheckMaker("satname");
  static oops::FilterMaker<OBS, oops::ObsFilter<OBS, ufo::ObsBoundsCheck> >
           boundsCheckMaker("Bounds Check");
  static oops::FilterMaker<OBS, oops::ObsFilter<OBS, ufo::BlackList> >
           blackListMaker("BlackList");
  static oops::FilterMaker<OBS, oops::ObsFilter<OBS, ufo::BlackList> >
           rejectListMaker("RejectList");  // alternative name
  static oops::FilterMaker<OBS, oops::ObsFilter<OBS, ufo::BackgroundCheck> >
           backgroundCheckMaker("Background Check");
  static oops::FilterMaker<OBS, oops::ObsFilter<OBS, ufo::BayesianBackgroundCheck> >
           BayesianBackgroundCheckMaker("Bayesian Background Check");
  static oops::FilterMaker<OBS, oops::ObsFilter<OBS, ufo::DifferenceCheck> >
           differenceCheckMaker("Difference Check");
  static oops::FilterMaker<OBS, oops::ObsFilter<OBS, ufo::HistoryCheck> >
           historyCheckMaker("History Check");
  static oops::FilterMaker<OBS, oops::ObsFilter<OBS, ufo::ModelObThreshold> >
           ModelObThresholdMaker("ModelOb Threshold");
  static oops::FilterMaker<OBS, oops::ObsFilter<OBS, ufo::ROobserror> >
           ROobserrorMaker("ROobserror");
  static oops::FilterMaker<OBS, oops::ObsFilter<OBS, ufo::Thinning> >
           thinningMaker("Thinning");
  static oops::FilterMaker<OBS, oops::ObsFilter<OBS, ufo::Gaussian_Thinning> >
           gaussianThinningMaker("Gaussian Thinning");
  static oops::FilterMaker<OBS, oops::ObsFilter<OBS, ufo::MWCLWCheck> >
           MWCLWCheckMaker("MWCLW Check");
  static oops::FilterMaker<OBS, oops::ObsFilter<OBS, ufo::ObsDomainErrCheck> >
           domainErrCheckMaker("DomainErr Check");
  static oops::FilterMaker<OBS, oops::ObsFilter<OBS, ufo::ProfileConsistencyChecks> >
           profileConsistencyChecksMaker("Profile Consistency Checks");
  static oops::FilterMaker<OBS, oops::ObsFilter<OBS, ufo::BackgroundCheckRONBAM> >
           backgroundCheckRONBAMMaker("Background Check RONBAM");
  static oops::FilterMaker<OBS, oops::ObsFilter<OBS, ufo::TemporalThinning> >
           temporalThinningMaker("Temporal Thinning");
  static oops::FilterMaker<OBS, oops::ObsFilter<OBS, ufo::PoissonDiskThinning> >
           poissonDiskThinningMaker("Poisson Disk Thinning");
  static oops::FilterMaker<OBS, oops::ObsFilter<OBS, ufo::ObsDiagnosticsWriter> >
           YDIAGsaverMaker("YDIAGsaver");
  static oops::FilterMaker<OBS, oops::ObsFilter<OBS, ufo::TrackCheck> >
           TrackCheckMaker("Track Check");
  static oops::FilterMaker<OBS, oops::ObsFilter<OBS, ufo::MetOfficeBuddyCheck> >
           MetOfficeBuddyCheckMaker("Met Office Buddy Check");
  static oops::FilterMaker<OBS, oops::ObsFilter<OBS, ufo::ObsDerivativeCheck> >
           DerivativeCheckMaker("Derivative Check");
  static oops::FilterMaker<OBS, oops::ObsFilter<OBS, ufo::TrackCheckShip> >
           ShipTrackCheckMaker("Ship Track Check");
  static oops::FilterMaker<OBS, oops::ObsFilter<OBS, ufo::StuckCheck> >
           StuckCheckMaker("Stuck Check");
  static oops::FilterMaker<OBS, oops::ObsFilter<OBS, ufo::GNSSROOneDVarCheck> >
           GNSSROOneDVarCheckMaker("GNSS-RO 1DVar Check");
  static oops::FilterMaker<OBS, oops::ObsFilter<OBS, ufo::VariableAssignment> >
           variableAssignmentMaker("Variable Assignment");
  static oops::FilterMaker<OBS, oops::ObsFilter<OBS, ufo::VariableTransforms> >
           VariableTransformsMaker("Variable Transforms");
  static oops::FilterMaker<OBS, oops::ObsFilter<OBS, ufo::ProfileBackgroundCheck> >
           ProfileBackgroundCheckMaker("Profile Background Check");
  static oops::FilterMaker<OBS, oops::ObsFilter<OBS, ufo::ProfileFewObsCheck> >
           ProfileFewObsCheckMaker("Profile Few Observations Check");
  static oops::FilterMaker<OBS, oops::ObsFilter<OBS, ufo::AcceptList> >
           acceptListMaker("AcceptList");
  static oops::FilterMaker<OBS, oops::ObsFilter<OBS, ufo::PerformAction> >
           performActionMaker("Perform Action");
  static oops::FilterMaker<OBS, oops::ObsFilter<OBS, ufo::BayesianBackgroundQCFlags> >
           BayesianBackgroundQCFlagsMaker("Bayesian Background QC Flags");
  static oops::FilterMaker<OBS, oops::ObsFilter<OBS, ufo::ProbabilityGrossErrorWholeReport> >
           ProbabilityGrossErrorWholeReportMaker("Bayesian Whole Report");
  static oops::FilterMaker<OBS, oops::ObsFilter<OBS, ufo::ImpactHeightCheck> >
           ImpactHeightCheckMaker("GNSSRO Impact Height Check");
  static oops::FilterMaker<OBS, oops::ObsFilter<OBS, ufo::SatwindInversionCorrection> >
             SatwindInversionCorrectionMaker("Satwind Inversion Correction");

  // Only include this filter if rttov is present
  #if defined(RTTOV_FOUND)
    static oops::FilterMaker<OBS, oops::ObsFilter<OBS, ufo::RTTOVOneDVarCheck> >
             RTTOVOneDVarCheckMaker("RTTOV OneDVar Check");
  #endif

  // For backward compatibility, register some filters under legacy names used in the past
  static oops::FilterMaker<OBS, oops::ObsFilter<OBS, ufo::Gaussian_Thinning> >
           legacyGaussianThinningMaker("Gaussian_Thinning");
  static oops::FilterMaker<OBS, oops::ObsFilter<OBS, ufo::TemporalThinning> >
           legacyTemporalThinningMaker("TemporalThinning");
}

}  // namespace ufo

#endif  // UFO_INSTANTIATEOBSFILTERFACTORY_H_
