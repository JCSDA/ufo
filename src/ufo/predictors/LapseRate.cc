/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <iterator>
#include <string>

#include "ufo/predictors/LapseRate.h"

#include "ioda/ObsSpace.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

static PredictorMaker<LapseRate> makerFuncLapseRate_("lapse_rate");

// -----------------------------------------------------------------------------

LapseRate::LapseRate(const eckit::Configuration & conf,
                     const std::vector<int> & jobs,
                     const std::string & sensor,
                     const eckit::mpi::Comm & comm)
  : PredictorBase(conf, jobs, sensor, comm), order_(1),
    tlaps_(), tsum_(), ntlapupdate_()
{
  // get the order if it is provided in options
  order_ = conf.getInt("predictor.options.order", 1);

    // override the predictor name for differentiable
  if (order_ > 1) {
    name() = name() + "_order_" + std::to_string(order_);
  }

  // required variables
  geovars_ += oops::Variables({"air_temperature",
                               "air_pressure",
                               "average_surface_temperature_within_field_of_view"});
  if (jobs_.size() > 0) {
    hdiags_ += oops::Variables({"transmittances_of_atmosphere_layer"}, jobs);
  } else {
    oops::Log::error() << "Channels size is ZERO !" << std::endl;
    ABORT("Channels size is ZERO !");
  }

  const auto tlapse_file = conf.getString("predictor.options.datain", "");

  const auto root = 0;

  if (!tlapse_file.empty()) {
    if (comm_.rank() == root) {
      ObsBiasIO< Record > tlapseIO(tlapse_file, std::ios::in);
      const auto tlapName = conf.getString("predictor.options.tlapse.name", "tlap");
      tlaps_ =
        tlapseIO.readByChannels(sensor_, jobs_, tlapName);

      tsum_ =
        tlapseIO.readByChannels(sensor_, jobs_, "tsum");

      ntlapupdate_ =
        tlapseIO.readByChannels(sensor_, jobs_, "ntlapupdate");
    }

    comm_.broadcast(tlaps_, root);
  } else {
    oops::Log::error() << "tlapse file is not provided !" << std::endl;
    ABORT("tlapse file is not provided !");
  }
}

// -----------------------------------------------------------------------------

void LapseRate::write(const eckit::Configuration & conf,
                      ObsBiasIO< Record > & fileOut) {
  const auto tlapName = conf.getString("predictor.options.tlapse.name", "tlap");
  fileOut.addByChannels(sensor_, jobs_, tlapName, tlaps_);
  fileOut.addByChannels(sensor_, jobs_, "tsum", tsum_);
  fileOut.addByChannels(sensor_, jobs_, "ntlapupdate", ntlapupdate_);
}

// -----------------------------------------------------------------------------

void LapseRate::compute(const ioda::ObsSpace & odb,
                        const GeoVaLs & geovals,
                        const ObsDiagnostics & ydiags,
                        ioda::ObsVector & out) const {
  const std::size_t njobs = jobs_.size();
  const std::size_t nlocs = odb.nlocs();

  // assure shape of out
  ASSERT(out.nlocs() == nlocs);

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
    hdiags = "transmittances_of_atmosphere_layer_" + std::to_string(jobs_[jb]);
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

  for (std::size_t jl = 0; jl < nlocs; ++jl) {
    for (std::size_t jb = 0; jb < njobs; ++jb) {
        tlapchn = (ptau5[jb][nlevs-2][jl]-ptau5[jb][nlevs-1][jl])*(tsavg5[jl]-tvp[nlevs-2][jl]);
        for (std::size_t k = 1; k < nlevs-1; ++k) {
          tlapchn = tlapchn+(ptau5[jb][nlevs-k-2][jl]-ptau5[jb][nlevs-k-1][jl])*
                    (tvp[nlevs-k][jl]-tvp[nlevs-k-2][jl]);
        }
        out[jl*njobs+jb] = pow((tlapchn - tlaps_[jb]), order_);
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
