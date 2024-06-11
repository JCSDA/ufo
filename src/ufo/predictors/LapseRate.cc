/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#include <fstream>
#include <string>
#include <vector>

#include "ufo/predictors/LapseRate.h"

#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

static PredictorMaker<LapseRate> makerFuncLapseRate_("lapseRate");

// -----------------------------------------------------------------------------

LapseRate::LapseRate(const Parameters_ & parameters, const oops::ObsVariables & vars)
  : PredictorBase(parameters, vars),
    order_(parameters.order.value().value_or(1))
{
  if (parameters.order.value() != boost::none) {
    // override the predictor name to distinguish between lapse_rate predictors of different orders
    name() = name() + "_order_" + std::to_string(order_);
  }

  // required variables
  geovars_ += oops::Variables({"air_temperature",
                               "air_pressure",
                               "average_surface_temperature_within_field_of_view"});
  if (vars.size() > 0) {
    hdiags_ += oops::ObsVariables({"transmittances_of_atmosphere_layer"}, vars.channels());
  } else {
    oops::Log::error() << "Channels size is ZERO !" << std::endl;
    ABORT("Channels size is ZERO !");
  }

  // This is a very preliminary method, please revisit
  // more flexibilites are needed
  const std::string & tlapse_file = parameters.tlapse;
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
    ABORT("Unable to open tlapse file ");
  }
}

// -----------------------------------------------------------------------------

void LapseRate::compute(const ioda::ObsSpace & odb,
                        const GeoVaLs & geovals,
                        const ObsDiagnostics & ydiags,
                        const ObsBias &,
                        ioda::ObsVector & out) const {
  const std::size_t nvars = out.nvars();
  const std::size_t nlocs = out.nlocs();

  // common vectors storage
  std::vector <float> pred(nlocs, 0.0);

  // retrieve the average surface temperature
  std::vector<float> tsavg5(nlocs, 0.0);
  geovals.get(tsavg5, oops::Variable{"average_surface_temperature_within_field_of_view"});

  // Retrieve the transmittances_of_atmosphere_layer from Hdiag
  std::vector<std::vector<std::vector<float>>> ptau5;
  std::vector<std::vector<float>> tmpvar;

  std::string hdiags;
  for (std::size_t jvar = 0; jvar < nvars; ++jvar) {
    hdiags = "transmittances_of_atmosphere_layer_" + std::to_string(vars_.channels()[jvar]);
    tmpvar.clear();
    for (std::size_t js = 0; js < ydiags.nlevs(hdiags); ++js) {
      ydiags.get(pred, hdiags, js);
      tmpvar.push_back(pred);
    }
    ptau5.push_back(tmpvar);
  }

  // Retrieve the temperature
  std::vector<std::vector<float>> tvp;
  std::size_t nlevs = geovals.nlevs(oops::Variable{"air_temperature"});
  for (std::size_t js = 0; js < nlevs; ++js) {
    geovals.getAtLevel(pred, oops::Variable{"air_temperature"}, js);
    tvp.push_back(pred);
  }
  nlevs = geovals.nlevs(oops::Variable{"air_pressure"});
  float tlapchn;

  // sort out the tlapmean based on vars
  std::vector<float> tlap;
  for (std::size_t jvar = 0; jvar < nvars; ++jvar) {
    auto it = tlapmean_.find(vars_.channels()[jvar]);
    if (it != tlapmean_.end()) {
      tlap.push_back(it->second);
    } else {
      oops::Log::error() << "Could not locate tlapemean for channel: " <<
                            vars_.channels()[jvar] << std::endl;
      ABORT("Could not locate tlapemean value");
    }
  }

  for (std::size_t jloc = 0; jloc < nlocs; ++jloc) {
    for (std::size_t jvar = 0; jvar < nvars; ++jvar) {
        tlapchn = (ptau5[jvar][nlevs-2][jloc]-ptau5[jvar][nlevs-1][jloc])*
                  (tsavg5[jloc]-tvp[nlevs-2][jloc]);
        for (std::size_t k = 1; k < nlevs-1; ++k) {
          tlapchn = tlapchn+(ptau5[jvar][nlevs-k-2][jloc]-ptau5[jvar][nlevs-k-1][jloc])*
                    (tvp[nlevs-k][jloc]-tvp[nlevs-k-2][jloc]);
        }
        out[jloc*nvars+jvar] = pow((tlapchn - tlap[jvar]), order_);
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
