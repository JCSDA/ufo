/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_WINDSSPDBCHECK_H_
#define UFO_FILTERS_OBSFUNCTIONS_WINDSSPDBCHECK_H_

#include <string>
#include <vector>

#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ufo {

///
///  This function require fives parameters: 'error_min' and
/// 'error_max' are lower and upper bound for 'ObsError'. 'cgross' is
///  the upper bound for gross-check, wndtype is the report type from the
///  PREPBUFR.  'variable' is assigned as either windEastward or windNorthward.
///  'ObsErrorData' (optional) is the specification of group name for
///  'ObsError',assigned to ObsErrorData by default.
///  In addition, there is an optional parameter HofX group name that can be
///  specified to override the default test_hofx group name used for testing purposes.
///  This can be set to GsiHofX, for example.
///
class WindsSPDBCheckParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(WindsSPDBCheckParameters, Parameters)

 public:
  /// the existence of min,max error values are required
  oops::RequiredParameter<std::vector<int>> wndtype{"wndtype", this};
  oops::RequiredParameter<std::vector<float>> error_min{"error_min", this};
  oops::RequiredParameter<std::vector<float>> error_max{"error_max", this};
  oops::RequiredParameter<std::vector<float>> cgross{"cgross", this};
  oops::RequiredParameter<std::string> checkvars{"variable", this};
  /// Name of the HofX group used to replace the default (HofX) group
  oops::Parameter<std::string> test_hofx{"test_hofx", "HofX", this};
  /// Name of the ObsError group used to replace the default (ObsErrorData) group
  oops::Parameter<std::string> testObserr{"test_obserr", "ObsErrorData", this};
  /// Name of the data group to which the QC flag is applied  (default is QCflagsData)
  oops::Parameter<std::string> testQCflag{"test_qcflag", "QCflagsData", this};
};

// -----------------------------------------------------------------------------

///  Compute the wind speed background difference (SPDB) of the observation and model.
///  If the SPDB is less than zero, then cgross = 0.7*cgross when the reported pressure
///  is bewteen 300hPa and 400hPa for PREPBUFR report type 244 and 245, and between
///  200hPa and 400hPa for PREPBUFR report type of 253 and 254. Next, compute a final
///  value of
///     cgross*ObsError*abs(OMF of eastward/northward wind)/(wind speed residual).
///  The final value should be compared against a function absolute threshold value in
///  background Check filter to flag obs locations as bad.  In regular usage, the
///  test_hofx option would be omitted in order for HofX values are used for model wind
///  components.  Also, an optional group name for the ObsError variable
///  can be supplied if different from ObsErrorData, and its value is bounded by
///  parameters error_min and error_max.
///
///
///        The PREPBUFR report type(wndtype) is described as folllows:

///        240  NESDIS IR(SHORT-WAVE) CLOUD DRIFT (ALL LEVELS) (GOES) - u, v
///        241  INDIA IR (LONG-WAVE) AND VISIBLE CLOUD DRIFT (ALL LEVELS)
///                                                  (INSAT, KALPANA) - u, v
///        242  JMA IR (LONG-WAVE) AND VISIBLE CLOUD DRIFT AT LEVELS
///                               BELOW 850 MB (GMS, MTSAT, HIMAWARI) - u, v
///        243  EUMETSAT IR (LONG-WAVE) AND VISIBLE CLOUD DRIFT AT LEVELS
///                                           BELOW 850 MB (METEOSAT) - u, v
///        244  AVHRR/POES IR (LONG-WAVE) CLOUD DRIFT (ALL LEVELS)
///                                                     (NOAA, METOP) - u, v
///        245  NESDIS IR (LONG-WAVE) CLOUD DRIFT (ALL LEVELS) (GOES) - u, v
///        246  NESDIS IMAGER WATER VAPOR (ALL LEVELS) - CLOUD TOP (GOES) - u, v
///        247  NESDIS IMAGER WATER VAPOR (ALL LEVELS) - DEEP LAYER (GOES) - u, v
///        248  NESDIS SOUNDER WATER VAPOR (ALL LEVELS) - CLOUD TOP (GOES) - u, v
///        249  NESDIS SOUNDER WATER VAPOR (ALL LEVELS) - DEEP LAYER (GOES) - u, v
///        250  JMA IMAGER WATER VAPOR (ALL LEVELS) - CLOUD TOP & DEEP LAYER
///                                                  (GMS, MTSAT, HIMAWARI) - u, v
///        251  NESDIS VISIBLE CLOUD DRIFT (ALL LEVELS) (GOES) - u, v
///        252  JMA IR (LONG-WAVE) AND VISIBLE CLOUD DRIFT AT LEVELS ABOVE 850 MB
///                                                  (GMS, MTSAT, HIMAWARI) - u, v
///        253  EUMETSAT IR (LONG-WAVE) AND VISIBLE CLOUD DRIFT AT LEVELS ABOVE
///                                                       850 MB (METEOSAT) - u, v
///        254  EUMETSAT IMAGER WATER VAPOR (ALL LEVELS) - CLOUD TOP & DEEP
///                                                        LAYER (METEOSAT) - u, v
///        255  NESDIS PICTURE TRIPLET CLOUD DRIFT (LOW LEVELS) (GOES)
///        256  INDIA IMAGER WATER VAPOR (ALL LEVELS) (INSAT, KALPANA) - u, v
///        257  MODIS/POES IR (LONG-WAVE) CLOUD DRIFT (ALL LEVELS)
///                                                       (AQUA, TERRA) - u,v
///        258  MODIS/POES IMAGER WATER VAPOR (ALL LEVELS) - CLOUD TOP
///                                                       (AQUA, TERRA) - u, v
///        259  MODIS/POES IMAGER WATER VAPOR (ALL LEVELS) - DEEP LAYER
///                                                       (AQUA, TERRA) - u, v
///        260  VIIRS/POES IR (LONG-WAVE) CLOUD DRIFT (ALL LEVELS) (NPP) - u,v
///
/// ### Sample YAML configuration for windEastward
/// - filter: Background Check
///   filter variables:
///   - name: windEastward
///    function absolute threshold:
///    - name: ObsFunction/WindsSPDBCheck
///      options:
///        wndtype: [  240,  241,  242,  243,  244,  245,  246,  247,  248,  249,
///                    250,  251,  252,  253,  254,  256,  257,  258,  259,  260]
///        cgross:  [  2.5,  2.5,  2.5,  1.5,  2.5,  1.3,  1.3,  2.5,  2.5,  2.5,
///                    2.5,  1.3,  2.5,  1.5,  1.5,  2.5,  2.5,  2.5,  2.5,  2.5]
///        error_min:  [1.4,  1.4,  1.4,  1.4,  1.4,  1.4,  1.4,  1.4,  1.4,  1.4,
///                     1.4,  1.4,  1.4,  1.4,  1.4,  1.4,  1.4,  1.4,  1.4,  1.4]
///        error_max:  [6.1,  6.1, 15.0, 15.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0,
///                    20.0, 20.0, 20.0, 20.0, 20.0, 20.1, 20.1, 20.1, 20.1, 20.1]
///        variable: windEastward
///    action:
//      name: reject
///
class WindsSPDBCheck : public ObsFunctionBase<float> {
 public:
  static const std::string classname() {return "WindsSPDBCheck";}

  explicit WindsSPDBCheck(const eckit::LocalConfiguration &
                              = eckit::LocalConfiguration());
  ~WindsSPDBCheck();

  void compute(const ObsFilterData &, ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  WindsSPDBCheckParameters options_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_WINDSSPDBCHECK_H_
