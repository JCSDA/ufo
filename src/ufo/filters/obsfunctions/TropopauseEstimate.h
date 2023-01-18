/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef UFO_FILTERS_OBSFUNCTIONS_TROPOPAUSEESTIMATE_H_
#define UFO_FILTERS_OBSFUNCTIONS_TROPOPAUSEESTIMATE_H_

#include <string>
#include <vector>

#include "oops/util/parameters/Parameter.h"

#include "ufo/filters/obsfunctions/ObsFunctionBase.h"
#include "ufo/filters/Variables.h"

namespace ioda {
  template<typename DATATYPE> class ObsDataVector;
}

namespace ufo {

class ObsFilterData;

///
/// \brief Two optional parameters permit overriding default tropopause pressure at the
///        equator and poles (linear interp between).  The optional save (default=false)
///        option can be used to save calculated tropopause estimate to the output file.
///        By default (convert_p2z=false), the output is tropopause pressure, but this
///        optional argument can convert the answer from pressure to height using the ICAO
///        standard atmosphere approximation.
///
class TropopauseEstimateParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(TropopauseEstimateParameters, Parameters)

 public:
  /// Default: The tropopause near the equator is located at a constant pressure level (135 hPa)
  /// in a belt from 15 S to 15 N and increases linearly to the poles to 360 hPa.
  oops::Parameter<float> tropo_equator{"tropo_equator", 13500.0f, this};
  oops::Parameter<float> tropo_pole{"tropo_pole", 36000.0f, this};
  oops::Parameter<bool> convert_p2z{"convert_p2z", false, this};
  oops::Parameter<bool> save{"save", false, this};
};

// -----------------------------------------------------------------------------

/// \brief Create a first-guess estimate of the tropopause pressure that is based
///        on latitude with some adjustment for day-of-year.  An optional parameter
///        can convert the final answer from pressure to height (convert_p2z: true).
///        The conversion is ultra simple approximation of ICAO std. atmosphere because
///        the code is making a tropopause *estimate* only.  Additional options can
///        alter the default tropopause pressure at equator (135 hPa) and poles (360 hPa).
///        The code in this method is crude and purely designed for estimating the tropopause
///        when lacking a model-derived estimate that may otherwise arrive via GeoVaLs.
///
/// ~~~
///
/// ### Sample YAML configuration
///     - filter: Difference Check
///       filter variables:
///       - name: windEastward
///       - name: windNorthward
///       reference: ObsFunction/TropopauseEstimate
///    #  options:                         # These options will not work yet with Difference Check
///    #    - tropo_equator: 13000         # 130 hPa
///    #    - tropo_pole: 37000            # 370 hPa
///       value: MetaData/pressure
///       minvalue: -5000                  # 50 hPa above tropopause level, negative p-diff
///
class TropopauseEstimate : public ObsFunctionBase<float> {
 public:
  static const std::string classname() {return "TropopauseEstimate";}

  explicit TropopauseEstimate(const eckit::LocalConfiguration &
                                 = eckit::LocalConfiguration());
  ~TropopauseEstimate();

  void compute(const ObsFilterData &, ioda::ObsDataVector<float> &) const;
  const ufo::Variables & requiredVariables() const;
 private:
  ufo::Variables invars_;
  TropopauseEstimateParameters options_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_FILTERS_OBSFUNCTIONS_TROPOPAUSEESTIMATE_H_
