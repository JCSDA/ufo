/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#include <fstream>
#include <iterator>
#include <map>
#include <string>
#include <vector>

#include "ufo/obsbias/predictors/LapseRate.h"

#include "ioda/ObsSpace.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

static PredictorMaker<LapseRate> makerFuncLapseRate_("lapse_rate");

// -----------------------------------------------------------------------------

LapseRate::LapseRate(const eckit::Configuration & conf)
  : PredictorBase(conf), order_(1)
{
  // get the order if it is provided in options
  if (conf.has("predictor.options.order")) {
    conf.get("predictor.options.order", order_);

    // override the predictor name for differentiable
    name() = name() + "_order_" + std::to_string(order_);
  }

  // required variables
  this->updateGeovars({"air_temperature",
                       "air_pressure",
                       "average_surface_temperature_within_field_of_view"});
  this->updateHdiagnostics({"transmittances_of_atmosphere_layer_CH"});

  // This is a very preliminary method, please revisit
  // more flexibilites are needed
  if (conf.has("predictor.options.tlapse")) {
    const std::string tlapse_file = conf.getString("predictor.options.tlapse");
    std::ifstream infile(tlapse_file);
    std::string nusis;   //  sensor/instrument/satellite
    int nuchan;  //  channel number
    float tlapse;

    if (infile.is_open()) {
      while (!infile.eof()) {
        infile >> nusis;
        infile >> nuchan;
        infile >> tlapse;
        tlapmean_[nuchan] = tlapse;
      }
      infile.close();
    } else {
      oops::Log::error() << "Unable to open file : "
                         << tlapse_file << std::endl;
      ABORT("Unable to open tlap file ");
    }
  } else {
    oops::Log::error() << "tlapse file is not provided !" << std::endl;
    ABORT("tlapse file is not provided !");
  }
}

// -----------------------------------------------------------------------------

void LapseRate::compute(const ioda::ObsSpace & odb,
                        const GeoVaLs & geovals,
                        const ObsDiagnostics & ydiags,
                        const std::vector<int> & jobs,
                        Eigen::MatrixXd & out) const {
  const std::size_t njobs = jobs.size();
  const std::size_t nlocs = odb.nlocs();

  // assure shape of out
  ASSERT(out.rows() == njobs && out.cols() == nlocs);

  // common vectors storage
  std::vector <float> pred(nlocs, 0.0);

  // retrieve the average surface temperature
  std::vector<float> tsavg5(nlocs, 0.0);
  geovals.get(tsavg5, "average_surface_temperature_within_field_of_view");

  // Retrieve the transmittances_of_atmosphere_layer from Hdiag
  std::vector<std::vector<std::vector<float>>> ptau5;
  std::vector<std::vector<float>> tmpvar;

  std::string hdiags;
  for (std::size_t jb = 0; jb < njobs; ++jb) {
    hdiags = "transmittances_of_atmosphere_layer_" + std::to_string(jobs[jb]);
    tmpvar.clear();
    for (std::size_t js = 0; js < ydiags.nlevs(hdiags); ++js) {
      ydiags.get(pred, hdiags, js+1);
      tmpvar.push_back(pred);
    }
    ptau5.push_back(tmpvar);
  }

  // Retrieve the temperature
  std::vector<std::vector<float>> tvp;
  std::size_t nlevs = geovals.nlevs("air_temperature");
  for (std::size_t js = 0; js < nlevs; ++js) {
    geovals.get(pred, "air_temperature", js+1);
    tvp.push_back(pred);
  }
  nlevs = geovals.nlevs("air_pressure");
  float tlapchn;

  // sort out the tlapmean based on jobs
  std::vector<float> tlap;
  for (std::size_t jb = 0; jb < njobs; ++jb) {
    auto it = tlapmean_.find(jobs[jb]);
    if (it != tlapmean_.end()) {
      tlap.push_back(it->second);
    } else {
      oops::Log::error() << "Could not locate tlapemean for channel: " << jobs[jb] << std::endl;
      ABORT("Could not locate tlapemean value");
    }
  }

  for (std::size_t jl = 0; jl < nlocs; ++jl) {
    for (std::size_t jb = 0; jb < njobs; ++jb) {
        tlapchn = (ptau5[jb][nlevs-2][jl]-ptau5[jb][nlevs-1][jl])*(tsavg5[jl]-tvp[nlevs-2][jl]);
        for (std::size_t k = 1; k < nlevs-1; ++k) {
          tlapchn = tlapchn+(ptau5[jb][nlevs-k-2][jl]-ptau5[jb][nlevs-k-1][jl])*
                    (tvp[nlevs-k][jl]-tvp[nlevs-k-2][jl]);
        }
        out(jb, jl) = pow((tlapchn - tlap[jb]), order_);
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
