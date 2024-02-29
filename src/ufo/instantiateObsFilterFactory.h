/*
 * (C) Copyright 2019-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_INSTANTIATEOBSFILTERFACTORY_H_
#define UFO_INSTANTIATEOBSFILTERFACTORY_H_

#include "oops/base/instantiateObsFilterFactory.h"
#include "oops/interface/ObsFilterBase.h"
#include "ufo/filters/AcceptList.h"
#include "ufo/filters/BackgroundCheck.h"
#include "ufo/filters/BayesianBackgroundCheck.h"
#include "ufo/filters/BayesianBackgroundQCFlags.h"
#include "ufo/filters/BlackList.h"
#include "ufo/filters/ConventionalProfileProcessing.h"
#include "ufo/filters/CopyFlagsFromExtendedToOriginalSpace.h"
#include "ufo/filters/CreateDiagnosticFlags.h"
#include "ufo/filters/DifferenceCheck.h"
#include "ufo/filters/FinalCheck.h"
#include "ufo/filters/Gaussian_Thinning.h"
#include "ufo/filters/gnssroonedvarcheck/GNSSROOneDVarCheck.h"
#include "ufo/filters/HistoryCheck.h"
#include "ufo/filters/ImpactHeightCheck.h"
#include "ufo/filters/MetOfficeBuddyCheck.h"
#include "ufo/filters/MetOfficeDuplicateCheck.h"
#include "ufo/filters/MetOfficePressureConsistencyCheck.h"
#include "ufo/filters/ModelBestFitPressure.h"
#include "ufo/filters/ModelObThreshold.h"
#include "ufo/filters/MWCLWCheck.h"
#include "ufo/filters/ObsBoundsCheck.h"
#include "ufo/filters/ObsDerivativeCheck.h"
#include "ufo/filters/ObsDiagnosticsWriter.h"
#include "ufo/filters/ObsDomainCheck.h"
#include "ufo/filters/ObsDomainErrCheck.h"
#include "ufo/filters/ObsRefractivityGradientCheck.h"
#include "ufo/filters/PerformAction.h"
#include "ufo/filters/PoissonDiskThinning.h"
#include "ufo/filters/PreQC.h"
#include "ufo/filters/PrintFilterData.h"
#include "ufo/filters/ProbabilityGrossErrorWholeReport.h"
#include "ufo/filters/ProcessAMVQI.h"
#include "ufo/filters/ProfileAverageObsToModLevels.h"
#include "ufo/filters/ProfileBackgroundCheck.h"
#include "ufo/filters/ProfileFewObsCheck.h"
#include "ufo/filters/ProfileMaxDifferenceCheck.h"
#include "ufo/filters/ProfileUnFlagObsCheck.h"
#include "ufo/filters/QCmanager.h"
#include "ufo/filters/SatName.h"
#include "ufo/filters/SatwindInversionCorrection.h"
#include "ufo/filters/SpikeAndStepCheck.h"
#include "ufo/filters/StuckCheck.h"
#include "ufo/filters/SuperOb.h"
#include "ufo/filters/SuperRefractionCheckImpactParameter.h"
#include "ufo/filters/SuperRefractionCheckNBAM.h"
#include "ufo/filters/TemporalThinning.h"
#include "ufo/filters/Thinning.h"
#include "ufo/filters/TrackCheck.h"
#include "ufo/filters/TrackCheckShip.h"
#include "ufo/filters/VariableAssignment.h"
#include "ufo/filters/VariableTransforms.h"
#include "ufo/operators/gnssro/QC/BackgroundCheckRONBAM.h"
#include "ufo/operators/gnssro/QC/ROobserror.h"

#if defined(GSW_FOUND)
  #include "ufo/filters/OceanVerticalStabilityCheck.h"
#endif

#if defined(RTTOV_FOUND)
  #include "ufo/filters/rttovonedvarcheck/RTTOVOneDVarCheck.h"
#endif

#include "ufo/ObsTraits.h"

namespace ufo {
void instantiateObsFilterFactory() {
  oops::instantiateObsFilterFactory<ObsTraits>();
  static oops::interface::FilterMaker<ObsTraits, AcceptList>
           acceptListMaker("AcceptList");
  static oops::interface::FilterMaker<ObsTraits, BackgroundCheck>
           backgroundCheckMaker("Background Check");
  static oops::interface::FilterMaker<ObsTraits, BackgroundCheckRONBAM>
           backgroundCheckRONBAMMaker("Background Check RONBAM");
  static oops::interface::FilterMaker<ObsTraits, BayesianBackgroundCheck>
           BayesianBackgroundCheckMaker("Bayesian Background Check");
  static oops::interface::FilterMaker<ObsTraits, BayesianBackgroundQCFlags>
           BayesianBackgroundQCFlagsMaker("Bayesian Background QC Flags");
  static oops::interface::FilterMaker<ObsTraits, ProbabilityGrossErrorWholeReport>
           ProbabilityGrossErrorWholeReportMaker("Bayesian Whole Report");
  static oops::interface::FilterMaker<ObsTraits, BlackList>
           blackListMaker("BlackList");  // same as RejectList
  static oops::interface::FilterMaker<ObsTraits, ObsBoundsCheck>
           boundsCheckMaker("Bounds Check");
  static oops::interface::FilterMaker<ObsTraits, ConventionalProfileProcessing>
           conventionalProfileProcessingMaker("Conventional Profile Processing");
  static oops::interface::FilterMaker<ObsTraits, CreateDiagnosticFlags>
           CreateDiagnosticFlagsMaker("Create Diagnostic Flags");
  static oops::interface::FilterMaker<ObsTraits, CopyFlagsFromExtendedToOriginalSpace>
           CopyFlagsFromExtendedToOriginalSpaceMaker("Copy Flags From Extended To Original Space");
  static oops::interface::FilterMaker<ObsTraits, ObsDerivativeCheck>
           DerivativeCheckMaker("Derivative Check");
  static oops::interface::FilterMaker<ObsTraits, DifferenceCheck>
           differenceCheckMaker("Difference Check");
  static oops::interface::FilterMaker<ObsTraits, ObsDomainCheck>
           domainCheckMaker("Domain Check");
  static oops::interface::FilterMaker<ObsTraits, ObsDomainErrCheck>
           domainErrCheckMaker("DomainErr Check");
  static oops::interface::FilterMaker<ObsTraits, FinalCheck>
           finalCheckMaker("Final Check");
  static oops::interface::FilterMaker<ObsTraits, Gaussian_Thinning>
           gaussianThinningMaker("Gaussian Thinning");
  static oops::interface::FilterMaker<ObsTraits, GNSSROOneDVarCheck>
           GNSSROOneDVarCheckMaker("GNSS-RO 1DVar Check");
  static oops::interface::FilterMaker<ObsTraits, ImpactHeightCheck>
           ImpactHeightCheckMaker("GNSSRO Impact Height Check");
  static oops::interface::FilterMaker<ObsTraits, HistoryCheck>
           historyCheckMaker("History Check");
  static oops::interface::FilterMaker<ObsTraits, MetOfficeBuddyCheck>
           MetOfficeBuddyCheckMaker("Met Office Buddy Check");
  static oops::interface::FilterMaker<ObsTraits, MetOfficeDuplicateCheck>
           MetOfficeDuplicateCheckMaker("Met Office Duplicate Check");
  static oops::interface::FilterMaker<ObsTraits, MetOfficePressureConsistencyCheck>
           MetOfficePressureConsistencyCheckMaker("Met Office Pressure Consistency Check");
  static oops::interface::FilterMaker<ObsTraits, ModelBestFitPressure>
           ModelBestFitPressureMaker("Model Best Fit Pressure");
  static oops::interface::FilterMaker<ObsTraits, ModelObThreshold>
           ModelObThresholdMaker("ModelOb Threshold");
  static oops::interface::FilterMaker<ObsTraits, MWCLWCheck>
           MWCLWCheckMaker("MWCLW Check");
  static oops::interface::FilterMaker<ObsTraits, ObsRefractivityGradientCheck>
           ObsRefractivityGradientCheckMaker("Obs Refractivity Gradient Check");
  static oops::interface::FilterMaker<ObsTraits, PerformAction>
           performActionMaker("Perform Action");
  static oops::interface::FilterMaker<ObsTraits, PoissonDiskThinning>
           poissonDiskThinningMaker("Poisson Disk Thinning");
  static oops::interface::FilterMaker<ObsTraits, PreQC>
           preQCMaker("PreQC");
  static oops::interface::FilterMaker<ObsTraits, PrintFilterData>
           printFilterDataMaker("Print Filter Data");
  static oops::interface::FilterMaker<ObsTraits, ProcessAMVQI>
             ProcessAMVQIMaker("Process AMV QI");
  static oops::interface::FilterMaker<ObsTraits, ProfileAverageObsToModLevels>
           ProfileAverageObsToModLevelsMaker("Average Observations To Model Levels");
  static oops::interface::FilterMaker<ObsTraits, ProfileBackgroundCheck>
           ProfileBackgroundCheckMaker("Profile Background Check");
  static oops::interface::FilterMaker<ObsTraits, ProfileFewObsCheck>
           ProfileFewObsCheckMaker("Profile Few Observations Check");
  static oops::interface::FilterMaker<ObsTraits, ProfileMaxDifferenceCheck>
           ProfileMaxDifferenceCheckMaker("Profile Max Difference Check");
  static oops::interface::FilterMaker<ObsTraits, ProfileUnFlagObsCheck>
           ProfileUnFlagObsCheckMaker("Profile Unflag Observations Check");
  static oops::interface::FilterMaker<ObsTraits, BlackList>
           rejectListMaker("RejectList");  // same as BlackList
  static oops::interface::FilterMaker<ObsTraits, ROobserror>
           ROobserrorMaker("ROobserror");
  static oops::interface::FilterMaker<ObsTraits, QCmanager>
           qcManagerMaker("QCmanager");
  static oops::interface::FilterMaker<ObsTraits, SatName>
           satnameCheckMaker("satname");
  static oops::interface::FilterMaker<ObsTraits, SatwindInversionCorrection>
             SatwindInversionCorrectionMaker("Satwind Inversion Correction");
  static oops::interface::FilterMaker<ObsTraits, TrackCheckShip>
           ShipTrackCheckMaker("Ship Track Check");
  static oops::interface::FilterMaker<ObsTraits, SpikeAndStepCheck>
           SpikeAndStepCheckMaker("Spike and Step Check");
  static oops::interface::FilterMaker<ObsTraits, StuckCheck>
           StuckCheckMaker("Stuck Check");
  static oops::interface::FilterMaker<ObsTraits, SuperRefractionCheckImpactParameter>
           SuperRefractionCheckImpactParameterCheckMaker("Impact Parameter Check");
  static oops::interface::FilterMaker<ObsTraits, SuperOb>
           SuperObMaker("SuperOb");
  static oops::interface::FilterMaker<ObsTraits, SuperRefractionCheckNBAM>
           SuperRefractionCheckNBAMCheckMaker("Super Refraction Check NBAM");
  static oops::interface::FilterMaker<ObsTraits, TemporalThinning>
           temporalThinningMaker("Temporal Thinning");
  static oops::interface::FilterMaker<ObsTraits, Thinning>
           thinningMaker("Thinning");
  static oops::interface::FilterMaker<ObsTraits, TrackCheck>
           TrackCheckMaker("Track Check");
  static oops::interface::FilterMaker<ObsTraits, VariableAssignment>
           variableAssignmentMaker("Variable Assignment");
  static oops::interface::FilterMaker<ObsTraits, VariableTransforms>
           VariableTransformsMaker("Variable Transforms");
  static oops::interface::FilterMaker<ObsTraits, ObsDiagnosticsWriter>
           YDIAGsaverMaker("YDIAGsaver");

  // Only include this filter if gsw is present
  #if defined(GSW_FOUND)
  static oops::interface::FilterMaker<ObsTraits, OceanVerticalStabilityCheck>
           OceanVerticalStabilityCheckMaker("Ocean Vertical Stability Check");
  #endif

  // Only include this filter if rttov is present
  #if defined(RTTOV_FOUND)
    static oops::interface::FilterMaker<ObsTraits, RTTOVOneDVarCheck>
             RTTOVOneDVarCheckMaker("RTTOV OneDVar Check");
  #endif

  // For backward compatibility, register some filters under legacy names used in the past
  static oops::interface::FilterMaker<ObsTraits, Gaussian_Thinning>
           legacyGaussianThinningMaker("Gaussian_Thinning");
  static oops::interface::FilterMaker<ObsTraits, TemporalThinning>
           legacyTemporalThinningMaker("TemporalThinning");
}

}  // namespace ufo

#endif  // UFO_INSTANTIATEOBSFILTERFACTORY_H_
