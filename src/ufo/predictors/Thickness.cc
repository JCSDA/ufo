/*
 * (C) Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>
#include <vector>

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "ufo/GeoVaLs.h"
#include "ufo/predictors/Thickness.h"
#include "ufo/utils/Constants.h"

namespace ufo {

static PredictorMaker<Thickness> makerFuncThickness_("thickness");

// -----------------------------------------------------------------------------

Thickness::Thickness(const Parameters_ & parameters, const oops::ObsVariables & vars)
    : PredictorBase(parameters, vars) {
  // required variables
  geovars_ += oops::Variables({oops::Variable{"air_temperature"},
                               oops::Variable{"air_pressure"}});
  // required options
  parameters_ = parameters;

  // override the predictor name to distinguish between
  // thickness predictors at different pressure layers
  name() = name() + "_" +  std::to_string(static_cast<int>(parameters_.layerBase.value()/100))
                  + '_' +  std::to_string(static_cast<int>(parameters_.layerTop.value()/100))
                  + "hPa";
}

// -----------------------------------------------------------------------------

void Thickness::compute(const ioda::ObsSpace & odb,
                       const GeoVaLs & geovals,
                       const ObsDiagnostics &,
                       const ObsBias &,
                       ioda::ObsVector & out) const {
  const std::size_t nlocs = odb.nlocs();
  const std::size_t nvars = out.nvars();

  // assure shape of out
  ASSERT(out.nlocs() == nlocs);

  const int t_levs = geovals.nlevs(oops::Variable{"air_temperature"});
  const int p_levs = geovals.nlevs(oops::Variable{"air_pressure"});
  std::vector <double> p_prof(t_levs, 0.0);
  std::vector<double> t_prof(p_levs, 0.0);
  std::vector <double> thick(odb.nlocs(), 0.0);

  const double p_high = parameters_.layerTop.value();
  const double p_low = parameters_.layerBase.value();
  const double pred_mean = parameters_.mean.value();
  const double pred_std_inv = 1.0/parameters_.stDev.value();;

  for (std::size_t jl = 0; jl < nlocs; ++jl) {
    geovals.getAtLocation(p_prof, oops::Variable{"air_pressure"}, jl);
    geovals.getAtLocation(t_prof, oops::Variable{"air_temperature"}, jl);

    if (p_prof.front() > p_prof.back())
      throw eckit::BadValue("GeoVaLs are in the wrong order.", Here());

    // Check that layer top is within pressure levels
    if (p_high > p_prof.back()) {
      oops::Log::error() << "layer top is greater than largest model pressure level" << std::endl;
      throw eckit::BadValue("layer top is greater than largest model pressure level", Here());
    }

    auto lower = std::lower_bound(p_prof.begin(), p_prof.end(), p_high);
    int i = lower - p_prof.begin();

    // find the thickness of the top fraction layer.
    double dp = p_prof[i] - p_high;
    double f = dp/(p_prof[i] - p_prof[i-1]);
    double p_av = p_prof[i] - 0.5*dp;
    double t_av = 0.5*((2.0-f)*t_prof[i]+ f*t_prof[i-1]);
    thick[jl] = t_av*(dp/p_av);
    i += 1;

    // For all the levels in the desired thickness band
    while (p_prof[i] < p_low && i < p_levs - 1) {
      dp = p_prof[i] - p_prof[i-1];
      p_av = p_prof[i] - 0.5*dp;
      t_av = 0.5*(t_prof[i] + t_prof[i-1]);
      thick[jl] += t_av*(dp/p_av);
      i += 1;
    }

    // When pressure level below thickness band or at last level
    if (p_prof[i] >= p_low) {
      dp = p_low - p_prof[i-1];
      f = dp/(p_prof[i] - p_prof[i-1]);
      p_av = p_low - 0.5*dp;
      t_av = 0.5*((f)*t_prof[i]+ (2.0-f)*t_prof[i-1]);
      thick[jl] += t_prof[i-1]*(dp/p_av);
      // Note that the line replicates OPS which contains a bug
      // the eqaution should be thick[jl] += t_av*(dp/p_av);
    } else {
      // Add thickness to last model level
      dp = p_prof[i] - p_prof[i-1];
      p_av = p_prof[i] - 0.5*dp;
      t_av = 0.5*(t_prof[i] + t_prof[i-1]);
      thick[jl] += t_av*(dp/p_av);
      // Then assume constant temperature until p_low
      dp = p_low - p_prof[i];
      p_av = p_low - 0.5*dp;
      t_av = t_prof[i];
      thick[jl] += t_av*(dp/p_av);
    }
  }

  const double km_per_m       = 1e-3;
  const double dry_air_gas_const = 287.0;  // Constants::rd not used for compatibility with OPS

  for (std::size_t jl = 0; jl < nlocs; ++jl) {
    for (std::size_t jb = 0; jb < nvars; ++jb) {
      out[jl*nvars+jb] = (thick[jl]
                         *(dry_air_gas_const/Constants::grav)*km_per_m - pred_mean)
                         *pred_std_inv;
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
