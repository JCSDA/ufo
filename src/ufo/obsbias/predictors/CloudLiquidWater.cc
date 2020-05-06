/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/obsbias/predictors/CloudLiquidWater.h"

#include <string>
#include <vector>

#include "ioda/ObsSpace.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

#include "ufo/filters/obsfunctions/CLWRetMW.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/utils/Constants.h"

namespace ufo {

static PredictorMaker<CloudLiquidWater>
       makerFuncCloudLiquidWater_("cloud_liquid_water");

// -----------------------------------------------------------------------------

CloudLiquidWater::CloudLiquidWater(const eckit::Configuration & conf)
  : PredictorBase(conf)
{
  // required variables
  this->updateGeovars({"water_area_fraction", "average_surface_temperature_within_field_of_view"});

  // Get channels for CLW retrieval from options
  if (conf.has("predictor.options")) {
    const eckit::LocalConfiguration options(conf.getSubConfiguration("predictor.options"));
    CLWRetMWParameters spec_opts;
    spec_opts.deserialize(options);

    ch238 = spec_opts.ch238.value();
    ch314 = spec_opts.ch314.value();
    ASSERT(spec_opts.ch238 !=0 && spec_opts.ch314 !=0);
  } else {
    oops::Log::error() << "CloudLiquidWater predictor needs options !" << std::endl;
    ABORT("CloudLiquidWater predictor needs options !");
  }
}

// -----------------------------------------------------------------------------

void CloudLiquidWater::compute(const ioda::ObsSpace & odb,
                               const GeoVaLs & geovals,
                               const ObsDiagnostics & ydiags,
                               const std::vector<int> & jobs,
                               Eigen::MatrixXd & out) const {
  const std::size_t njobs = jobs.size();
  const std::size_t nlocs = odb.nlocs();

  // assure shape of out
  ASSERT(out.rows() == njobs && out.cols() == nlocs);

  std::vector<int> channels(jobs.begin(), jobs.end());

  // Retrieve the brightness temperature from ODB
  std::vector<float> bt238(nlocs), bt314(nlocs);
  odb.get_db("ObsValue", "brightness_temperature_" + std::to_string(ch238), bt238);
  odb.get_db("ObsValue", "brightness_temperature_" + std::to_string(ch314), bt314);

  // Retrieve the water fraction from geovals
  std::vector<float> water_frac(nlocs, 0.0);
  geovals.get(water_frac, "water_area_fraction");

  std::vector<float> clw(nlocs);

  // retrieve the average surface temperature
  std::vector<float> tsavg(nlocs, 0.0);
  geovals.get(tsavg, "average_surface_temperature_within_field_of_view");

  // retreive the zenith angle
  std::vector<float> szas(nlocs, 0.0);
  odb.get_db("MetaData", "sensor_zenith_angle", szas);

  // compute the cloud liquid water
  CLWRetMW::cloudLiquidWater(szas, tsavg, water_frac, bt238, bt314, clw, nlocs);

  // weighted by cos(zenith_angle)
  for (std::size_t jl = 0; jl < nlocs; ++jl) {
    out.col(jl).setConstant(clw[jl]*cos(szas[jl])*cos(szas[jl]));
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
