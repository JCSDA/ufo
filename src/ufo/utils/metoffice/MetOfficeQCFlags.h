/*
 * (C) Crown copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_UTILS_METOFFICE_METOFFICEQCFLAGS_H_
#define UFO_UTILS_METOFFICE_METOFFICEQCFLAGS_H_

namespace ufo {
  /// \file QC flags used in OPS
  /// A variety of flags are defined for entire observations,
  /// particular elements within (generic) observations, and specific observations

namespace MetOfficeQCFlags {
  // Report flags for whole observations
  // Notes:
  //  PermRejectReport   = used for station list rejections
  //  SurplusReport      = used for ship/airep data judged to be (almost) duplicate
  //                       used for thinned buoy reports
  //  OutOfAreaReport    = used for data outside model area
  //                       (outside ocean mask for ocean/sst/wave analysis)
  //                       used for data outside time window of analysis
  enum WholeObReport {
    FinalRejectReport     = 1 << 0,  ///< One of flags 1-6 set
    PermRejectReport      = 1 << 1,  ///< Blacklisted data
    SurplusReport         = 1 << 2,  ///< (Near) duplicate data
    OutOfAreaReport       = 1 << 3,  ///< Outside analysis area/time
    LandRejectReport      = 1 << 4,  ///< Marine ob over land
    UsedInSuperObReport   = 1 << 5,  ///< Combined with other obs
    TrackRejectReport     = 1 << 6,  ///< Failed track check
    SuperObReport         = 1 << 8,  ///< This is a SuperOb
    // Sonde data
    NoPressureSensor      = 1 << 13,  ///< No PILOT pressure sensor
    // Satellite wind data
    MissingDataReport     = 1 << 14,  ///< Missing data
    SatwindAltReport      = 1 << 15,  ///< Satwind alternative p/uv
    SatwindGoodConstraint = 1 << 16,  ///< Best-fit pressure is not well constrained
    // Other miscellaneous flags
    Thin4DFlag            = 1 << 17,  ///< Duplicate found
    StationListThinFlag   = 1 << 18   ///< Rejected based on station list
  };

  // Assim flags for whole observations
  enum WholeObAssim {
    NewReport          = 1 << 0,  ///< Ob not yet assimilated
    AssimilatedReport  = 1 << 1,  ///< Ob already assimilated
  };

  // Flags for individual (generic) observation elements
  enum Elem {
    NoAssimFlag        = 1 << 23,  ///< Do not use in analysis
    FinalRejectFlag    = 1 << 0,   ///< Final QC flag
    BuddyRejectFlag    = 1 << 1,   ///< PGE>0.5 after buddy check
    BackRejectFlag     = 1 << 2,   ///< PGE>0.5 after backgr check
    PermRejectFlag     = 1 << 3,   ///< Blacklisted data
    ClimRejectFlag     = 1 << 4,   ///< PGE>0.5 after climat check
    BuddyPerfFlag      = 1 << 5,   ///< Buddy check performed
    BackPerfFlag       = 1 << 6,   ///< Background check performed
    ClimPerfFlag       = 1 << 7,   ///< Climatological check perf
    PermCorrectFlag    = 1 << 8,   ///< Fixed correction
    DataCorrectFlag    = 1 << 9,   ///< Eg sign correction
    ConsistencyFlag    = 1 << 10,  ///< Internal consistency check
    ExtremeValueFlag   = 1 << 11   ///< Extreme value check
  };

  // Flags for surface data
  enum Surface {
    TendencyFlag       = 1 << 12,  ///< Pressure tendency check.
    PstdRepFlag        = 1 << 13,  ///< Pstd reported not Pmsl.
    PstnPrefFlag       = 1 << 14,  ///< Use Pstn if reported.
    PmslUsedFlag       = 1 << 15,  ///< Pmsl used in P* calc.
    PstdUsedFlag       = 1 << 16,  ///< Pstd used in P* calc.
    PstnUsedFlag       = 1 << 17,  ///< Pstn used in P* calc.
    QNHinHgFlag        = 1 << 16,  ///< QNH in 0.01 inches Hg
    QNHhPaFlag         = 1 << 17,  ///< QNH in whole hPa - Metars
    RHreportFlag       = 1 << 18,  ///< RH was reported
    SiteQualityFlag    = 1 << 20,  ///< AIRQAL site quality reject flag
    VisRejFlag         = 1 << 22,  ///< Reject Visibility Ob
    notRoundedFlag     = 1 << 24,  ///< Metar QNH not rounded to whole hPa
  };

  // Flags for profiles
  enum Profile {
    HydrostaticFlag    = 1 << 12,  ///< Hydrostatic check flag
    InterpolationFlag  = 1 << 13,  ///< Interpolation check flag
    SuperadiabatFlag   = 1 << 14,  ///< Superadiabatic check flag
    SurfaceLevelFlag   = 1 << 15,  ///< Surface Level
    StandardLevelFlag  = 1 << 16,  ///< Standard Level
    SigTempLevelFlag   = 1 << 17,  ///< Significant Temperature
    SigWindLevelFlag   = 1 << 18,  ///< Significant Wind Level
    MaxWindLevelFlag   = 1 << 19,  ///< Maximum Wind Level
    TropopauseFlag     = 1 << 20,  ///< Tropopause Level
    PartialLayerFlag   = 1 << 21   ///< Partial Layer Vert Average
  };

  // Flags for satellite winds
  enum SatWind {
    SatwindConfFlag       = 1 << 12,  ///< Satwind product confidence
    SatwindInversionFlag  = 1 << 13,  ///< Inversion height corrected
    SatwindDryLayerFlag   = 1 << 14,  ///< Model dry layer QC
    SatwindWrongLayerFlag = 1 << 15   ///< Wrong moist layer QC
  };

  // Flags for scatterometers
  enum Scatt {
    ScatConfidenceFlag  = 1 << 12,  ///< Wind product confidence
    ScatAmbigRemov1Flag = 1 << 13,  ///< Wind ambiguity removal
    ScatAmbigRemov2Flag = 1 << 14,  ///< Wind ambiguity removal
    ScatIncAngle1Flag   = 1 << 15,  ///< Wind angle of incidence
    ScatIncAngle2Flag   = 1 << 16   ///< Wind angle of incidence
  };

  // Flags for aircraft relative humidity
  enum AircraftRH {
    DerivedFromMixRatioFlag    = 1 << 12,  ///< Relative humidity derived from mixing ratio
    DerivedFromFlightLevelFlag = 1 << 13,  ///< Pressure derived from flight level
  };

  // Flags for satellite SST
  enum SatSST {
    DaytimeFlag     = 1 << 12,  ///< Observation recorded in daytime
    DiurnalWarmFlag = 1 << 13   ///< Indicates a likely diurnal warming component in signal
  };

  // Profile flags which depend on sounding data type
  enum Sounding {
    // TEMP and PILOT
    TEMPSigWind    = 1 << 1,   ///< Significant wind level
    TEMPSigTemp    = 1 << 2,   ///< Significant temperature level
    TEMPMaxWind    = 1 << 3,   ///< Maximum wind level
    TEMPTropopause = 1 << 4,   ///< Tropopause level
    TEMPStandard   = 1 << 5,   ///< Standard level
    TEMPSurface    = 1 << 6,   ///< Surface level
    TEMPStandardX  = 1 << 7,   ///< Semi-standard level
    // BUFR
    BUFRSigWind    = 1 << 11,  ///< Significant wind level
    BUFRSigTemp    = 1 << 13,  ///< Significant temperature level
    BUFRMaxWind    = 1 << 14,  ///< Maximum wind level
    BUFRTropopause = 1 << 15,  ///< Tropopause level
    BUFRStandard   = 1 << 16,  ///< Standard level
    BUFRSurface    = 1 << 17,  ///< Surface level
    BUFRStandardX  = 1 << 16   ///< Semi-standard level, grouped with standard levels in this case
  };
}  // namespace MetOfficeQCFlags

}  // namespace ufo

#endif  // UFO_UTILS_METOFFICE_METOFFICEQCFLAGS_H_
