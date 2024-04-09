/*
 * (C) Copyright 2017-2023 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef UFO_OPERATORS_RADARREFLECTIVITY_DIRECTZDA_OBSDIRECTZDA_H_
#define UFO_OPERATORS_RADARREFLECTIVITY_DIRECTZDA_OBSDIRECTZDA_H_

#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"

#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"

#include "ufo/ObsOperatorBase.h"
#include "ufo/ObsOperatorParametersBase.h"
#include "ufo/operators/radarreflectivity/directZDA/ObsDirectZDA.interface.h"
#include "ufo/utils/parameters/ParameterTraitsVariable.h"

/// Forward declarations
namespace ioda {
  class ObsSpace;
  class ObsVector;
}

namespace ufo {
  class GeoVaLs;
  class ObsDiagnostics;

  enum class MicrophysicsOption {
    THOMPSON, NSSL, LIN, GFDL
  };

  struct MicrophysicsOptionParameterTraitsHelper {
    typedef MicrophysicsOption EnumType;
    static constexpr char enumTypeName[] = "MicrophysicsOption";
    static constexpr util::NamedEnumerator<MicrophysicsOption> namedValues[] = {
      { MicrophysicsOption::THOMPSON, "Thompson" },
      { MicrophysicsOption::NSSL, "NSSL"},
      { MicrophysicsOption::LIN, "Lin"},
      { MicrophysicsOption::GFDL, "GFDL"}
    };
  };
}  // namespace ufo

namespace oops {

template <>
struct ParameterTraits<ufo::MicrophysicsOption> :
    public EnumParameterTraits<ufo::MicrophysicsOptionParameterTraitsHelper>
{};

}  // namespace oops

namespace ufo {


/// Configuration options recognized by the radar-reflectivity directZDA operator.

class ObsDirectZDAParameters: public ObsOperatorParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsDirectZDAParameters, ObsOperatorParametersBase)

 public:
  oops::Parameter<std::string> VertCoord
    {"vertical coordinate",
     "vertical coordinate used by the model",
     "geopotential_height",  // this should be consistent with var_prs defined in ufo_vars_mod
     this};

  oops::Parameter<MicrophysicsOption> microOption
    {"microphysics option",
     "microphysics option by name (e.g. Thompson, NSSL, GFDL, others)",
     MicrophysicsOption::THOMPSON,
     this};

  oops::Parameter<bool> use_variational
    {"use variational method",
     "use variational TL/AD method (direct ZDA)",
     false,
     this};

  /*
   * The following list of hydrometeor species mixing ratios and number concentrations
   * need to be consistent with list of variables in ufo_variables_mod, which has some
   * confusing names that need to be resolved in the future - all to match CCPP
   * convention.  In the meantime, FV3 and MPAS use different variables.  The code
   * will default to FV3 names (rain_water, snow_water, and graupel), but all input
   * species can be overridden with these optional YAML parameters, which for MPAS
   * is mass_content_of_rain_in_atmosphere_layer for example.
   */

  oops::Parameter<std::string> var_rain_mixing_ratio
    {"var_rain_mixing_ratio",
     "Name of model rain mixing ratio variable",
     "rain_water",  // this should be consistent with var_qr
     this};

  oops::Parameter<std::string> var_snow_mixing_ratio
    {"var_snow_mixing_ratio",
     "Name of model snow mixing ratio variable",
     "snow_water",  // this should be consistent with var_qs
     this};

  oops::Parameter<std::string> var_graupel_mixing_ratio
    {"var_graupel_mixing_ratio",
     "Name of model graupel mixing ratio variable",
     "graupel",  // this should be consistent with var_qg
     this};

  oops::Parameter<std::string> var_hail_mixing_ratio
    {"var_hail_mixing_ratio",
     "Name of model hail mixing ratio variable",
     "hail",  // this should be consistent with var_qh
     this};

  oops::Parameter<std::string> var_rain_number_concentration
    {"var_rain_number_concentration",
     "Name of model rain number concentration variable",
     "rain_number_concentration",  // this should be consistent with var_nr
     this};

  oops::Parameter<std::string> var_snow_number_concentration
    {"var_snow_number_concentration",
     "Name of model snow number concentration variable",
     "snow_number_concentration",  // this should be consistent with var_ns
     this};

  oops::Parameter<std::string> var_graupel_number_concentration
    {"var_graupel_number_concentration",
     "Name of model graupel number concentration variable",
     "graupel_number_concentration",  // this should be consistent with var_ng
     this};

  oops::Parameter<std::string> var_hail_number_concentration
    {"var_hail_number_concentration",
     "Name of model hail number concentration variable",
     "hail_number_concentration",  // this should be consistent with var_nh
     this};

  oops::Parameter<std::string> var_graupel_vol_mixing_ratio
    {"var_graupel_vol_mixing_ratio",
     "Name of model graupel volume mixing ratio variable",
     "volume_mixing_ratio_of_graupel_in_air",  // this should be consistent with var_qvg
     this};

  oops::Parameter<std::string> var_hail_vol_mixing_ratio
    {"var_hail_vol_mixing_ratio",
     "Name of model hail volume mixing ratio variable",
     "volume_mixing_ratio_of_hail_in_air",  // this should be consistent with var_qvh
     this};
};

/*
 *   This radar reflectivity code has two internal components, one is the older more
 * traditional reflectivity calculated using the integral of [m(D)]^2 (mass-squared)
 * times N(D) (number distribution assumption) per unit diameter, D; where m(D) is
 * the assumed mass as function of diameter.  The second method is referred to as
 * directZDA which allows the hydrometeor mixing ratios (and number concentrations)
 * to be scaled and used with Tangent Linear and Adjoint models as described in the
 * following three references:
 *   Liu, C., M. Xue, and R. Kong, 2020: Direct variational assimilation of radar
 * reflectivity and radial velocity data: Issues with nonlinear reflectivity operator
 * and solutions. Mon. Wea Rev., 148, 1483-1502, https://doi.org/10.1175/MWR-D-19-0149.1.
 *   Chen, L., C. Liu, M. Xue, R. Kong, and Y. Jung, 2021: Use of power transform mixing
 * ratios as hydrometeor control variables for direct assimilation of radar reflectivity
 * in GSI En3DVar and tests with five convective storms cases. Mon. Wea. Rev., 149,
 * 645-659. https://doi.org/10.1175/MWR-D-20-0149.1.
 *   Li, H., C. Liu, M. Xue, J. Park, L. Chen, Y. Jung, R. Kong, and C.-C. Tong, 2022: Use
 * of power transform total number concentration as control variable for direct
 * assimilation of radar reflectivity in GSI En3DVar and tests with six convective storms
 * cases. Mon. Wea. Rev., 150, 821-842. https://doi.org/10.1175/MWR-D-21-0041.1.
 */

// -----------------------------------------------------------------------------
/// RadarReflectivity observation operator class
class ObsDirectZDA: public ObsOperatorBase,
                   private util::ObjectCounter<ObsDirectZDA> {
 public:
  typedef ioda::ObsDataVector<int> QCFlags_t;
  typedef ObsDirectZDAParameters Parameters_;
  static const std::string classname() {return "ufo::ObsDirectZDA";}

  ObsDirectZDA(const ioda::ObsSpace &, const Parameters_ &);
  virtual ~ObsDirectZDA();

// Obs Operator
  void simulateObs(const GeoVaLs &, ioda::ObsVector &, ObsDiagnostics &,
                   const QCFlags_t &) const override;

// Other
  const oops::Variables & requiredVars() const override {return varin_;}

  int & toFortran() {return keyOperDirectZDA_;}
  const int & toFortran() const {return keyOperDirectZDA_;}

 private:
  void print(std::ostream &) const override;
  F90hop keyOperDirectZDA_;
  const ioda::ObsSpace& odb_;
  oops::Variables varin_;
};

// -----------------------------------------------------------------------------

}  // namespace ufo

#endif  // UFO_OPERATORS_RADARREFLECTIVITY_DIRECTZDA_OBSDIRECTZDA_H_
