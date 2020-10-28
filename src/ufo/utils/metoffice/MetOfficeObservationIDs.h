/*
 * (C) Crown copyright 2020, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_UTILS_METOFFICE_METOFFICEOBSERVATIONIDS_H_
#define UFO_UTILS_METOFFICE_METOFFICEOBSERVATIONIDS_H_

namespace ufo {
  /// \file Numerical observation identifiers used in OPS.
  /// The observations are grouped into enums according to their category.

namespace MetOfficeObsIDs {

    // Surface observations
    enum Surface {
      Synop           = 10100,  ///< generic synop
      SynopManual     = 10101,  ///< non-automatic synop
      SynopAuto       = 10102,  ///< automatic synop
      Ship            = 10200,  ///< generic ship
      ShipManual      = 10201,  ///< non-automatic ship
      ShipAuto        = 10202,  ///< automatic ship
      PlatsRiggs      = 10204,  ///< platform/rigg
      ShipMrdBuoy     = 10210,  ///< moored buoy in ship code
      ShpMrdBySpOb    = 10211,  ///< moored buoy super ob from ship code reports
      ShipDrifter     = 10212,  ///< drifting buoy in ship code
      ShpDrftSpOb     = 10213,  ///< drifting buoy super ob from ship code reports
      Buoy            = 10300,  ///< generic buoy
      MooredBuoy      = 10310,  ///< moored buoy in buoy code
      MrdBuoySpOb     = 10311,  ///< moored buoy super ob
      Drifter         = 10312,  ///< drifting buoy in buoy code
      DrfBuoySpOb     = 10313,  ///< drifting buoy super ob
      Precip          = 10400,  ///< Radar precipitation
      Gauge           = 10401,  ///< Rain gauge accumulation
      Srew            = 10500,  ///< SREW precipitation
      RadWind         = 10601,  ///< Radar radial winds
      RadRefl         = 10602,  ///< Radar reflectivity
      RadTRH          = 10603,  ///< Radar 1DVar
      RadopWind       = 10604,  ///< Radar radial winds
      Radop2Wind      = 10605,  ///< Radar radial winds (RADOP2)
      RadRefr         = 10606,  ///< Radar radial refractivity
      RadVel          = 10607,  ///< Radar radial winds
      Wavenet         = 10700,  ///< Wavenet ocean wave spectra
      SynopMob        = 10800,  ///< mobile synop
      OpenRoad        = 10900,  ///< OPENROAD data
      AirQal          = 11000,  ///< Air quality ob
      AirQalEU        = 11001,  ///< Air quality EU ob
      Metar           = 11100,  ///< generic metar
      MetarManual     = 11101,  ///< non-automatic metar
      MetarAuto       = 11102,  ///< automatic metar
      Sferics         = 11200,  ///< SFERICS lightning strikes
      ATDNET          = 11300,  ///< SFERICS lightning strikes
      WOW             = 11400,  ///< WOW observations
      GroundLidar     = 11501,  ///< Lidarnet ceilometer backscatter
      SynopBufr       = 11600,  ///< generic synop
      BuoyBufr        = 11700,  ///< buoy bufr
      MooredBuoyBufr  = 11701,  ///< moored buoy bufr
      DrifterBuoyBufr = 11702   ///< drifting buoy bufr
    };

    // Satellite observations
    enum Satellite {
      Esaura         = 20100,  ///< ERS UARS(wave height)
      EsauraSpOb     = 20101,  ///< ERS UARS super ob
      Esauwa         = 20200,  ///< ERS ESA spectral wave
      SeaWinds       = 20300,  ///< Ku Band  scat winds
      AATSRUnkn      = 20400,  ///< ERS AATSR SST (unknown)
      AATSRDay       = 20401,  ///< ERS AATSR SST (day)
      AATSRNight     = 20402,  ///< ERS AATSR SST (night)
      Esauwi         = 20500,  ///< ERS scat winds
      UKMOSSTu       = 20600,  ///< UKMOSST unknown
      UKMOSSTd       = 20601,  ///< UKMOSST daytime
      UKMOSSTn       = 20602,  ///< UKMOSST nighttime
      UKMOSSTt       = 20603,  ///< UKMOSST twilight
      ALADIN         = 20700,  ///< HLOS winds from Aeolus(ALADIN)
      Satob          = 21000,  ///< Satob
      // Satob SSTs
      NoaaPo         = 21001,  ///< NOAA plar orbitor
      MeteoSat       = 21002,  ///< MeteoSat SST satob
      Himawari       = 21003,  ///< Himawari SST satob
      IndSat         = 21004,  ///< Indian SST satob
      GoesE          = 21005,  ///< GOES E SST satob
      GoesW          = 21006,  ///< GOES W SST satob
      AuoSat2        = 21007,  ///< AutoSat 2 SST
      // Satwinds
      ObsTYPEAMV     = 21010,  ///< Combined AMV data
      SatobIR        = 21021,  ///< Satob - Infrared
      SatobVIS       = 21022,  ///< Satob - Visible
      SatobWV        = 21023,  ///< Satob - Water Vapour
      SatobMC        = 21024,  ///< Satob - Multi-channel
      Goesamw        = 22000,  ///< HD GOES Satob, infrared
      Goesvis        = 22002,  ///< HD GOES Satob, visible
      Goeswv         = 22003,  ///< HD GOES Satob, water vapour
      Esahrvw        = 22100,  ///< EUMETSAT - visible
      Esacmw         = 22200,  ///< EUMETSAT - combined (LR)
      Esacmwir       = 22201,  ///< EUMETSAT - infrared (LR)
      Esacmwvis      = 22202,  ///< EUMETSAT - visible (LR)
      Esacmwwv       = 22203,  ///< EUMETSAT - water vapour (LR)
      Esacswvw       = 22300,  ///< EUMETSAT - clear sky water vapour
      Esahrwvw       = 22400,  ///< EUMETSAT - cloudy water vapour
      Goesbufr       = 22500,  ///< NESDIS - unknown channel
      GoesbufrIR     = 22501,  ///< NESDIS - infrared 10.8
      GoesbufrVIS    = 22502,  ///< NESDIS - visible
      GoesbufrWV     = 22503,  ///< NESDIS - cloudy water vapour
      GoesbufrCSWV   = 22505,  ///< NESDIS - clear sky water vapour
      GoesbufrWV62   = 22531,  ///< NESDIS - cloudy water vapour 6.2
      GoesbufrWV73   = 22532,  ///< NESDIS - cloudy water vapour 7.3
      GoesbufrCSWV62 = 22551,  ///< NESDIS - clear sky water vapour 6.2
      GoesbufrCSWV73 = 22552,  ///< NESDIS - clear sky water vapour 7.3
      GoesbufrIR38   = 22511,  ///< NESDIS - infrared 3.8
      GoesbufrIR16   = 22515,  ///< NESDIS - infrared 1.6
      Modis          = 22600,  ///< MODIS - unknown channel
      ModisIR        = 22601,  ///< MODIS - infrared 10.8
      ModisWV        = 22603,  ///< MODIS - cloudy water vapour
      ModisMIXWV     = 22604,  ///< MODIS - mixed water vapour
      ModisCSWV      = 22605,  ///< MODIS - clear sky water vapour
      ModisIR38      = 22611,  ///< MODIS - infrared 3.8
      ModisIR16      = 22615,  ///< MODIS - infrared 1.6
      Jmawinds       = 23500,  ///< JMA - unknown channel
      JmawindsIR     = 23501,  ///< JMA - infrared
      JmawindsVIS    = 23502,  ///< JMA - visible
      JmawindsWV     = 23503,  ///< JMA - cloudy water vapour
      JmawindsCSWV   = 23505,  ///< JMA - clear sky water vapour
      JmawindsMIXWV  = 23507,  ///< JMA - mixed water vapour
      JmawindsIR38   = 23511,  ///< JMA - infrared 3.8
      JmawindsWV62   = 23531,  ///< JMA - cloudy water vapour 6.2
      JmawindsWV73   = 23532,  ///< JMA - cloudy water vapour 7.3
      JmawindsCSWV62 = 23551,  ///< JMA - clear sky water vapour 6.2
      JmawindsCSWV73 = 23552,  ///< JMA - clear sky water vapour 7.3
      Msgwinds       = 23600,  ///< MSG - unknown channel
      MsgIR38        = 23611,  ///< MSG - infrared 3.8
      MsgIR87        = 23612,  ///< MSG - infrared 8.7
      MsgIR108       = 23613,  ///< MSG - infrared 10.8
      MsgOZONE       = 23614,  ///< MSG - infrared 9.7
      MsgHRVIS       = 23621,  ///< MSG - high res. visible
      MsgVIS06       = 23622,  ///< MSG - visible 0.6
      MsgVIS08       = 23623,  ///< MSG - visible 0.8
      MsgWV62        = 23631,  ///< MSG - cloudy water vapour 6.2
      MsgWV73        = 23632,  ///< MSG - cloudy water vapour 7.3
      MsgCSWV62      = 23651,  ///< MSG - clear sky water vapour 6.2
      MsgCSWV73      = 23652,  ///< MSG - clear sky water vapour 7.3
      Cmawinds       = 24900,  ///< CMA - unknown channel
      CmawindsIR     = 24901,  ///< CMA - infrared
      CmawindsVIS    = 24902,  ///< CMA - visible
      CmawindsWV     = 24903,  ///< CMA - cloudy water vapour
      CmawindsCSWV   = 24905,  ///< CMA - clear sky water vapour
      CmawindsMIXWV  = 24907,  ///< CMA - mixed water vapour
      CmawindsIR38   = 24911,  ///< CMA - infrared 3.8
      Imdwinds       = 25200,  ///< IMD - unknown channel
      ImdwindsIR     = 25201,  ///< IMD - infrared
      ImdwindsVIS    = 25202,  ///< IMD - visible
      ImdwindsWV     = 25203,  ///< IMD - cloudy WV
      ImdwindsCSWV   = 25205,  ///< IMD - clear sky WV
      ImdwindsMIXWV  = 25207,  ///< IMD - mixed WV
      ImdwindsIR38   = 25211,  ///< IMD - infrared 3.8
      Kmawinds       = 25900,  ///< KMA - unknown channel
      KmawindsIR     = 25901,  ///< KMA - infrared
      KmawindsVIS    = 25902,  ///< KMA - visible
      KmawindsWV     = 25903,  ///< KMA - cloudy WV
      KmawindsCSWV   = 25905,  ///< KMA - clear sky WV
      KmawindsMIXWV  = 25907,  ///< KMA - mixed WV
      KmawindsIR38   = 25911,  ///< KMA - infrared 3.8
      Ukwinds        = 26400,  ///< UK MSG - unknown channel
      UkIR38         = 26411,  ///< UK MSG - infrared 3.8
      UkIR87         = 26412,  ///< UK MSG - infrared 8.7
      UkIR108        = 26413,  ///< UK MSG - infrared 10.8
      UkIR120        = 26415,  ///< UK MSG - infrared 12.0
      UkOZONE        = 26414,  ///< UK MSG - infrared 9.7
      UkHRVIS        = 26421,  ///< UK MSG - high res. visible
      UkVIS06        = 26422,  ///< UK MSG - visible 0.6
      UkVIS08        = 26423,  ///< UK MSG - visible 0.8
      UkWV62         = 26431,  ///< UK MSG - cloudy water vapour 6.2
      UkWV73         = 26432,  ///< UK MSG - cloudy water vapour 7.3
      UkCSWV62       = 26451,  ///< UK MSG - clear sky water vapour 6.2
      UkCSWV73       = 26452,  ///< UK MSG - clear sky water vapour 7.3
      Stereomv       = 27100,  ///< STEREOMV - unknown type
      StereomvVIS06  = 27122,  ///< STEREOMV - MISR winds (VIS 0.6)
      TOVS_G         = 21100,  ///< global TOVS
      TOVS_G2        = 21101,  ///< products only
      TOVS_G1D       = 21102,  ///< radiances/products
      TOVS_G1E       = 21103,  ///< G1D with NESDIS rad's
      TOVS_L         = 21200,  ///< local TOVS
      TOVS_L2        = 21201,  ///< products only
      TOVS_L1D       = 21202,  ///< radiances/products
      ATOVS_G        = 21400,  ///< global ATOVS
      ATOVS_G2       = 21401,  ///< products only
      ATOVS_G1D      = 21402,  ///< radiances/products
      ATOVS_L        = 21500,  ///< local ATOVS
      ATOVS_L2       = 21501,  ///< products only
      ATOVS_L1D      = 21502,  ///< radiances/products
      MWSFY3B        = 21450,  ///< Global FY3B as ATOVS subtype
      Gpsiwv         = 21700,  ///< Ground GPS
      ATMS           = 21800,  ///< ATMS radiances
      ATMSHR         = 21900,  ///< High res ATMS radiances
      AIRS           = 22800,  ///< AIRS Radiances
      GPSRO          = 22900,  ///< GPSRO
      SAltSSH        = 23000,  ///< Altimeter SSH
      SAltSSHT       = 23001,  ///< Altimeter SSH (real time)
      SAltSSHS       = 23002,  ///< Altimeter SSH with extra variables
      SSMIS          = 23100,  ///< SSMI/S
      AMSUB          = 23200,  ///< full res AMSUB
      SBUVozone      = 23700,  ///< SBUV ozone
      AIRSRR         = 23900,  ///< AIRS Reconstructed Radiances
      AIRSWF         = 24000,  ///< AIRS Warmest FOV Radiances
      IASI_G         = 24100,  ///< Global IASI
      IASI_L         = 24200,  ///< Local  IASI
      ASCAT          = 24300,  ///< ASCAT BUFR
      GHRSST         = 24400,  ///< GHRSST
      SEVIRIEUM      = 24500,  ///< SEVIRI EUMETSAT radiance
      SEVIRIAUTO     = 24600,  ///< SEVIRI AUTOSAT radiance
      SEVIRIAUTOUK   = 24610,  ///< SEVIRI AUTOSAT radiance for UK area
      SEVIRIAUTOFD   = 24650,  ///< SEVIRI AUTOSAT radiance for full disc
      SEVIRICTP      = 24660,  ///< SEVIRI AUTOSAT MSGCTP data
      SEVIRICTPUK    = 24670,  ///< SEVIRI AUTOSAT MSGCTP data for UK area
      SEVIRICTPFD    = 24680,  ///< SEVIRI AUTOSAT MSGCTP data for full disc
      AIRSMS         = 24700,  ///< AIRS Clearest MODIS FOV
      WINDSAT        = 24800,  ///< WINDSAT BUFR
      AMSR2          = 25000,  ///< AMSR2 radiances
      AIRSL          = 25100,  ///< AIRS locally received data
      IASI_H         = 25300,  ///< High Res IASI
      MERIS          = 25400,  ///< MERIS
      OLCI           = 25410,  ///< OLCI
      ASCATHR        = 25500,  ///< HiRes ASCAT BUFR
      MSGAOD         = 25600,  ///< SEVIRI AOD
      MVIRICSR       = 25700,  ///< MVIRI radiances
      GOESImCSR      = 25800,  ///< GOES Imager radiances
      CrIS           = 26000,  ///< CrIS+ATMS (global, thinned)
      CrISFSR        = 26001,  ///< CrIS FSR +ATMS (global, thinned)
      CRISHR         = 26100,  ///< CrIS+ATMS (UK, full res)
      CRISFSRHR      = 26101,  ///< CrIS FSR +ATMS (UK, full res)
      COMSCSR        = 26200,  ///< COMS Imager radiances
      MTSATImCSR     = 26300,  ///< MTSAT Imager radiances
      AHICSR         = 26310,  ///< AHI radiances
      AHIASR         = 26320,  ///< AHI all-sky radiances
      ABICSR         = 26330,  ///< ABI (GOES-16 onwards) radiances
      SATAOD         = 26600,  ///< Satellite AOD (Polar)
      ASCATCO        = 26700,  ///< ASCAT coastal winds
      IN3DImCSR      = 26800,  ///< INSAT3D Imager radiances
      IN3DSdCSR      = 26900,  ///< INSAT3D Sounder radiances
      MWSFY3         = 27000,  ///< FY-3C/D MWTS-2 and MWHS-2 radiances
      SAPHIR         = 28000,  ///< MT SAPHIR radiances
      MWRI           = 29000,  ///< MWRI radiances
      GMIlow         = 29100,  ///< GMI low freq channels
      GMIhigh        = 29200   ///< GMI high freq channels
    };

    // Aircraft observations
    enum Aircraft {
      Amdar  = 30100,  ///< Amdar
      Airep  = 30200,  ///< airep
      Tamdar = 30300,  ///< Tamdar
      WISDOM = 30400,  ///< WISDOM
      ModeS  = 30500   ///< MODE-S
    };

    // Bogus observations
    enum Bogus {
      TCBogus = 40100,  ///< TC Bogus
      Bogus   = 40300   ///< Surface or sonde bogus
    };

    // Atmospheric profile observations
    enum AtmosphericProfile {
      Temp        = 50100,  ///< generic temp
      TempLand    = 50101,  ///< temp
      TempShip    = 50102,  ///< temp ship
      TempMobile  = 50103,  ///< temp mobile
      Pilot       = 50200,  ///< generic pilot
      PilotLand   = 50201,  ///< pilot
      PilotShip   = 50202,  ///< pilot ship
      PilotMobile = 50203,  ///< pilot mobile
      DropSonde   = 50300,  ///< drop sonde
      WindProf    = 50400,  ///< wind profiler
      Sonde       = 50500,  ///< sonde (BUFR)
      TSTSonde    = 50501   ///< sonde (BUFR)
    };

    // Ocean profile observations
    enum OceanProfile {
      Bathy    = 60100,  ///< Bathy
      Tesac    = 60200,  ///< Tesac
      Argo     = 60201,  ///< ARGO floats
      ArgoBufr = 60202,  ///< ARGO floats (BUFR)
      BuoyProf = 60300,  ///< Buoy profiles
      OceanFB  = 60404,  ///< NetCDF Ocean Ferrybox
      OceanRE  = 60402,  ///< NetCDF Ocean EXRE0206 EXRE0186
      OceanTS  = 60403   ///< NetCDF Ocean FNCM FAC8862 FNHO FHZI HOWN VHW5167
    };

    // Sea ice observations
    enum SeaIce {
      SeaIce   = 60500,  ///< Sea ice concentration (OSSEAICE)
      HRSeaIce = 60600,  ///< Sea ice concentration (HRSEAICE)
      SeaIceN  = 60700   ///< Sea ice concentration (OSSEAICN)
    };

}  // namespace MetOfficeObsIDs

}  // namespace ufo

#endif  // UFO_UTILS_METOFFICE_METOFFICEOBSERVATIONIDS_H_
