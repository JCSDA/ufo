/*
 * (C) Copyright 2024-2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 * Author: Rong Kong (NCAR/MMM) - C++ rewrite of the PPRO operator from Fortran.
 */

#include "ufo/operators/radarpolarimetric/ppro/ObsPPRO.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <ostream>
#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/GeoVaLs.h"
#include "ufo/operators/radarpolarimetric/ppro/TCWA2Forward.h"
#include "ufo/utils/Constants.h"
#include "ufo/utils/VertInterp.interface.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsPPRO> makerPPRO_("PPRO");
// -----------------------------------------------------------------------------

namespace {

// Operator output variable names, in the canonical order zh, zdr, kdp, phv.
const std::array<std::string, 4> kOpVarNames = {
    "equivalentReflectivityFactor",   // zh
    "differentialReflectivity",       // zdr
    "specificDifferentialPhase",      // kdp
    "copolarCorrelationCoefficient"   // phv
};

// Vertical-coordinate variable names (mirrors ufo_variables_mod).
const char *const kVarGeopotentialHeight = "geopotential_height";
const char *const kVarHeightAboveMSL = "height_above_mean_sea_level";

// Human-readable microphysics name for logging.
std::string microName(MicrophysicsOption opt) {
  switch (opt) {
    case MicrophysicsOption::THOMPSON: return "Thompson";
    case MicrophysicsOption::WSM6:     return "WSM6";
    case MicrophysicsOption::NSSL:     return "NSSL";
    case MicrophysicsOption::LIN:      return "Lin";
    case MicrophysicsOption::GFDL:     return "GFDL";
    case MicrophysicsOption::TCWA2:    return "TCWA2";
    default:                           return "unknown";
  }
}

}  // namespace

// -----------------------------------------------------------------------------

ppro::PPROTuning ObsPPRO::buildTuning(const ObsPPROParameters &params) {
  ppro::PPROTuning t;
  t.skip_small_qx = params.skip_small_qx.value();

  t.enable_melting_transition = params.enable_melting_transition.value();
  t.snow_ratio_low = params.snow_ratio_low.value();
  t.snow_ratio_high = params.snow_ratio_high.value();
  t.graupel_ratio_low = params.graupel_ratio_low.value();
  t.graupel_ratio_high = params.graupel_ratio_high.value();
  t.hail_ratio_low = params.hail_ratio_low.value();
  t.hail_ratio_high = params.hail_ratio_high.value();

  t.melting_water_fraction = params.melting_water_fraction.value();

  t.enable_dm_melting_limit = params.enable_dm_melting_limit.value();
  t.dm_melting_transition_width = params.dm_melting_transition_width.value();

  t.enable_dm_regularization = params.enable_dm_regularization.value();

  t.tuning_dm_rain = params.tuning_dm_rain.value();
  t.tuning_dm_melting_snow = params.tuning_dm_melting_snow.value();
  t.tuning_dm_melting_graupel = params.tuning_dm_melting_graupel.value();
  t.tuning_dm_melting_hail = params.tuning_dm_melting_hail.value();
  t.tuning_dm_pure_snow = params.tuning_dm_pure_snow.value();
  t.tuning_dm_pure_graupel = params.tuning_dm_pure_graupel.value();
  t.tuning_dm_pure_hail = params.tuning_dm_pure_hail.value();

  t.tuning_melt_frac_snow = params.tuning_melt_frac_snow.value();
  t.tuning_melt_frac_graupel = params.tuning_melt_frac_graupel.value();
  t.tuning_melt_frac_hail = params.tuning_melt_frac_hail.value();

  t.melting_scheme_option = params.melting_scheme.value();

  t.dmmax_rain = params.dmmax_rain.value();
  t.dmmax_pure_snow = params.dmmax_pure_snow.value();
  t.dmmax_melting_snow = params.dmmax_melting_snow.value();
  t.dmmax_pure_graupel = params.dmmax_pure_graupel.value();
  t.dmmax_melting_graupel = params.dmmax_melting_graupel.value();
  t.dmmax_pure_hail = params.dmmax_pure_hail.value();
  t.dmmax_melting_hail = params.dmmax_melting_hail.value();

  t.treat_hail_as_graupel = params.treat_hail_as_graupel.value();

  t.melting_exp_rain = params.melting_exp_rain.value();
  t.melting_exp_ice = params.melting_exp_ice.value();
  t.fraction_exp_rain = params.fraction_exp_rain.value();
  t.fraction_exp_ice = params.fraction_exp_ice.value();

  return t;
}

// -----------------------------------------------------------------------------

ObsPPRO::ObsPPRO(const ioda::ObsSpace & odb, const ObsPPROParameters & params)
  : ObsOperatorBase(odb), odb_(odb), varin_(),
    microOption_(params.microOption.value()),
    useTcwa2_(params.pproOption.value() == PolarimetricOperatorOption::TCWA2),
    coeffMelt_(params.coeff_melt.value()),
    vCoord_(params.VertCoord.value()),
    outColumn_{{-1, -1, -1, -1}},
    zhang21_(buildTuning(params))
{
  // Validate the configuration up front, before building the GeoVaLs list.
  if (useTcwa2_ && microOption_ != MicrophysicsOption::TCWA2) {
    throw eckit::BadValue(
        "ObsPPRO: TCWA2 polarimetric operator requires TCWA2 microphysics", Here());
  }

  // Validate the vertical coordinate.
  if (vCoord_ != kVarGeopotentialHeight && vCoord_ != kVarHeightAboveMSL &&
      vCoord_ != "height") {
    throw eckit::BadValue(
        "ObsPPRO: unsupported vertical coordinate '" + vCoord_ +
        "' (use geopotential_height, height_above_mean_sea_level, or height)",
        Here());
  }

  // Validate the melting scheme.  Only 'liu24' and 'jung08' are implemented;
  // reject anything else explicitly rather than silently falling back to liu24.
  const std::string & meltScheme = params.melting_scheme.value();
  if (meltScheme != "liu24" && meltScheme != "jung08") {
    throw eckit::BadValue(
        "ObsPPRO: unsupported melting scheme '" + meltScheme +
        "' (supported: liu24, jung08)", Here());
  }

  // Build the ordered list of required model variables (mirrors geovars_list in
  // ufo_PPRO_mod.F90).  The first four are fixed; the rest depend on the
  // microphysics scheme and the selected polarimetric operator.
  geovars_.push_back("moist_air_density");
  geovars_.push_back("air_temperature");
  geovars_.push_back("air_pressure");
  geovars_.push_back("water_vapor_mixing_ratio_wrt_moist_air");
  geovars_.push_back(params.var_rain_mixing_ratio.value());
  geovars_.push_back(params.var_snow_mixing_ratio.value());
  geovars_.push_back(params.var_graupel_mixing_ratio.value());

  switch (microOption_) {
    case MicrophysicsOption::THOMPSON:
      geovars_.push_back(params.var_rain_number_concentration.value());
      break;
    case MicrophysicsOption::WSM6:
      // 1-moment: rain/snow/graupel mixing ratios only, using the WSM6
      // single-moment PSD assumptions (assumed N0 intercepts) in dmZwsm6.
      break;
    case MicrophysicsOption::NSSL:
      geovars_.push_back(params.var_hail_mixing_ratio.value());
      geovars_.push_back(params.var_rain_number_concentration.value());
      geovars_.push_back(params.var_snow_number_concentration.value());
      geovars_.push_back(params.var_graupel_number_concentration.value());
      geovars_.push_back(params.var_hail_number_concentration.value());
      geovars_.push_back(params.var_graupel_vol_mixing_ratio.value());
      geovars_.push_back(params.var_hail_vol_mixing_ratio.value());
      break;
    case MicrophysicsOption::TCWA2:
      geovars_.push_back(params.var_rain_number_concentration.value());
      geovars_.push_back(params.var_snow_number_concentration.value());
      geovars_.push_back(params.var_graupel_number_concentration.value());
      if (useTcwa2_) {
        geovars_.push_back(params.var_cloud_mixing_ratio.value());
        geovars_.push_back(params.var_ice_mixing_ratio.value());
        geovars_.push_back(params.var_ice_number_concentration.value());
        geovars_.push_back(params.var_snow_melted_fraction.value());
        geovars_.push_back(params.var_graupel_melted_fraction.value());
      }
      break;
    default:
      throw eckit::BadValue("ObsPPRO: unsupported microphysics option", Here());
  }

  // Register required GeoVaLs: all model variables plus the vertical coordinate.
  for (const std::string &name : geovars_) {
    varin_.push_back(name);
  }
  varin_.push_back(vCoord_);

  // Map each simulated variable to one of the operator outputs, and record its
  // column index in the ObsVector.  Every simulated variable must be a known
  // PPRO output (mirrors the abort in ufo_PPRO_setup).
  const std::vector<std::string> & simVars = odb.assimvariables().variables();
  for (std::size_t i = 0; i < simVars.size(); ++i) {
    const auto it = std::find(kOpVarNames.begin(), kOpVarNames.end(), simVars[i]);
    if (it == kOpVarNames.end()) {
      throw eckit::BadValue(
          "ObsPPRO: simulated variable '" + simVars[i] +
          "' is not a PPRO operator variable", Here());
    }
    outColumn_[static_cast<std::size_t>(it - kOpVarNames.begin())] =
        static_cast<int>(i);
  }

  oops::Log::trace() << "ObsPPRO created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsPPRO::~ObsPPRO() {
  oops::Log::trace() << "ObsPPRO destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsPPRO::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                          ObsDiagnostics &,
                          const QCFlags_t & qc_flags) const {
  oops::Log::trace() << "ObsPPRO: observation operator starting" << std::endl;

  const std::size_t nlocs = odb_.nlocs();
  const int nvars = ovec.nvars();
  const double missing = util::missingValue<double>();

  const oops::Variable vcoordVar{vCoord_};
  const std::size_t nlevs = gv.nlevs(vcoordVar);

  // Observation vertical coordinate (always MetaData/height, as in Fortran).
  std::vector<double> obsHeight(nlocs);
  odb_.get_db("MetaData", "height", obsHeight);

  // Radar band: 1 = S-band, 2 = C-band; default to S-band if not provided.
  std::vector<int> iband(nlocs, 1);
  if (odb_.has("MetaData", "iband")) {
    odb_.get_db("MetaData", "iband", iband);
  }

  // Pre-build oops::Variable handles for each required model variable.
  std::vector<oops::Variable> geovarHandles;
  geovarHandles.reserve(geovars_.size());
  for (const std::string &name : geovars_) geovarHandles.emplace_back(name);

  const int nlev = static_cast<int>(nlevs);
  std::vector<double> vcoordProfile(nlevs);
  std::vector<double> profile(nlevs);
  std::vector<double> fields(geovars_.size());

  const double rd = Constants::rd;
  const double D608 = 0.608;

  auto setMissing = [&](std::size_t jloc) {
    for (int k = 0; k < 4; ++k) {
      if (outColumn_[k] >= 0) ovec[jloc * nvars + outColumn_[k]] = missing;
    }
  };

  for (std::size_t jloc = 0; jloc < nlocs; ++jloc) {
    const double obl = obsHeight[jloc];

    // Vertical interpolation weights from the model vertical coordinate, using
    // the shared UFO vert_interp utility (identical to the original operator).
    gv.getAtLocation(vcoordProfile, vcoordVar, jloc);
    int wi;
    double wf;
    vert_interp_weights_f90(nlev, obl, vcoordProfile.data(), wi, wf);

    // Interpolate each required model variable to the observation level.
    for (std::size_t v = 0; v < geovars_.size(); ++v) {
      gv.getAtLocation(profile, geovarHandles[v], jloc);
      vert_interp_apply_f90(nlev, profile.data(), fields[v], wi, wf);
    }

    const double t = fields[1];     // temperature [K]
    const double p = fields[2];     // pressure [Pa]
    double q = fields[3];           // specific humidity
    q = q / (1.0 - q);              // -> water vapor mixing ratio
    const double rho = p / (rd * t * (1.0 + D608 * q));  // dry air density

    if (rho < 0.0) {
      setMissing(jloc);
      continue;
    }

    const double qr = 1000.0 * fields[4];  // kg/kg -> g/kg
    const double qs = 1000.0 * fields[5];
    const double qg = 1000.0 * fields[6];

    double zh = 0.0, zdr = 1.0, kdp = 0.0, phv = 1.0;

    if (microOption_ == MicrophysicsOption::TCWA2 && useTcwa2_) {
      const double nr = fields[7] * rho;
      const double ns = fields[8] * rho;
      const double ng = fields[9] * rho;
      const double qc = 1000.0 * fields[10];
      const double qi = 1000.0 * fields[11];
      const double ni = fields[12] * rho;
      const double smlf = fields[13];
      const double gmlf = fields[14];
      ppro::tcwa2ComputePoint(iband[jloc], rho, t, qr, qs, qg, qc, qi,
                              nr, ns, ng, ni, smlf, gmlf, zh, zdr, kdp, phv);
    } else {
      ppro::Zhang21Optionals opt;
      // The melting tuning coefficient applies to every microphysics option
      // handled by the Zhang21 operator, not just NSSL.
      opt.coeff_melt = coeffMelt_;
      switch (microOption_) {
        case MicrophysicsOption::THOMPSON:
          opt.nr = fields[7] * rho;
          break;
        case MicrophysicsOption::WSM6:
          break;
        case MicrophysicsOption::NSSL:
          opt.qh = 1000.0 * fields[7];
          opt.nr = fields[8] * rho;
          opt.ns = fields[9] * rho;
          opt.ng = fields[10] * rho;
          opt.nh = fields[11] * rho;
          // fields[12]/fields[13] are the graupel/hail volume mixing ratios.
          // They are read as required GeoVaLs and reserved for an optional
          // variable-density treatment (opt.vg/opt.vh), but are intentionally
          // left unset so Zhang21 uses fixed graupel/hail density by default.
          opt.temperature = t - 273.15;
          break;
        case MicrophysicsOption::TCWA2:
          // TCWA2 microphysics handled by the Zhang21 operator (2-moment, no hail).
          opt.nr = fields[7] * rho;
          opt.ns = fields[8] * rho;
          opt.ng = fields[9] * rho;
          break;
        default:
          break;
      }
      zhang21_.computePoint(iband[jloc], rho, t, qr, qs, qg,
                            zh, zdr, kdp, phv, opt);
    }

    // Convert to output units and apply caps (mirrors ufo_PPRO_simobs).
    // Floor negative/zero reflectivity at 0 dBZ: weak echo (0 < zh < 1) yields
    // a negative dBZ, and no echo (zh <= 0) cannot be converted; both are set
    // to 0 dBZ so clear-air / weak-echo information is retained.
    const double zh_dBZ = (zh > 0.0) ? std::max(0.0, 10.0 * std::log10(zh)) : 0.0;

    // ZDR < 1 (linear) is unphysical for rain and gives a negative dB; floor it
    // at 0 dB. Values >= 1 are capped at 6 dB.
    double zdr_dB;
    if (zdr >= 1.0) {
      zdr_dB = std::min(6.0, 10.0 * std::log10(zdr));
    } else {
      zdr_dB = 0.0;
    }

    if (kdp > 0.0) kdp = std::min(6.0, kdp);

    const double out[4] = {zh_dBZ, zdr_dB, kdp, phv};
    for (int k = 0; k < 4; ++k) {
      if (outColumn_[k] >= 0) ovec[jloc * nvars + outColumn_[k]] = out[k];
    }
  }

  oops::Log::trace() << "ObsPPRO: observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsPPRO::print(std::ostream & os) const {
  os << "ObsPPRO: parameterized polarimetric radar operator" << std::endl;
  os << "  polarimetric operator = " << (useTcwa2_ ? "TCWA2" : "Zhang21")
     << std::endl;
  os << "  microphysics option   = " << microName(microOption_) << std::endl;
  os << "  vertical coordinate   = " << vCoord_ << std::endl;
  os << "  required GeoVaLs (" << geovars_.size() << "):";
  for (const std::string &name : geovars_) os << " " << name;
  os << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
