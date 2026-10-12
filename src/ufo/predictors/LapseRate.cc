/*
 * (C) Copyright 2020-2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#include <cmath>
#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"

#include "ufo/predictors/LapseRate.h"

#include "ioda/distribution/Accumulator.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsBias.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

static PredictorMaker<LapseRate> makerFuncLapseRate_("lapseRate");

// Shared per-channel mean lapse rates for all LapseRate instances on one ObsSpace.
// One static map per simulation is keyed to unique ObsSpace objects. This is typical
// for most real applications, and ensures numerical consistency between NL and TL
// in oops TL/AD tests, where some drivers hold two ObsBias objects for the same ObsSpace.
struct MeansSlot {
  bool ready = false;
  std::vector<float> means;
};
static std::map<const ioda::ObsSpace *, MeansSlot> meansByObsSpace_;

// -----------------------------------------------------------------------------

LapseRate::LapseRate(const Parameters_ & parameters, const oops::ObsVariables & vars)
  : PredictorBase(parameters, vars),
    useOnlineMeans_(parameters.onlineMean.value()),
    order_(parameters.order.value().value_or(1))
{
  const bool hasTlapseFile = parameters.tlapse.value() != boost::none;
  if (hasTlapseFile == useOnlineMeans_) {
    throw eckit::UserError(
      "LapseRate predictors must specify exactly one of 'tlapse' (file path) or "
      "'online mean: true'",
      Here());
  }

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

  if (!useOnlineMeans_) {
    // This is a very preliminary method, please revisit
    // more flexibilites are needed
    const std::string & tlapse_file = parameters.tlapse.value().value();
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
  std::vector<float> pred(nlocs, 0.0f);
  // retrieve the average surface temperature
  std::vector<float> tsavg5(nlocs, 0.0f);
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
  const std::size_t nlevs = geovals.nlevs(oops::Variable{"air_temperature"});
  for (std::size_t js = 0; js < nlevs; ++js) {
    geovals.getAtLevel(pred, oops::Variable{"air_temperature"}, js);
    tvp.push_back(pred);
  }

  const float fmiss = util::missingValue<float>();
  const double dmiss = util::missingValue<double>();
  const auto isBad = [fmiss](float v) { return v == fmiss || std::isnan(v); };

  std::vector<std::vector<float>> tlapchn_all(nvars, std::vector<float>(nlocs, fmiss));
  if (nvars == 0 || ptau5.empty() || ptau5[0].size() < nlevs || tvp.size() < nlevs ||
        nlevs < 2) {
    oops::Log::error() << "No valid transmittance profiles or air temperature profiles to compute lapse rates" << std::endl;
    ABORT("No valid transmittance profiles or air temperature profiles to compute lapse rates");
  } else {
    float tlapchn;
    for (std::size_t jloc = 0; jloc < nlocs; ++jloc) {
      for (std::size_t jvar = 0; jvar < nvars; ++jvar) {
        // Check for missing or NaN values.
        bool missingInput = isBad(tsavg5[jloc]);
        for (std::size_t lev = 0; lev < nlevs && !missingInput; ++lev) {
          missingInput = isBad(tvp[lev][jloc]) ||
                         isBad(ptau5[jvar][lev][jloc]);
        }
        if (missingInput) {
          continue;
        }

        tlapchn = (ptau5[jvar][nlevs-2][jloc]-ptau5[jvar][nlevs-1][jloc])*
                  (tsavg5[jloc]-tvp[nlevs-2][jloc]);
        for (std::size_t k = 1; k < nlevs-1; ++k) {
          tlapchn = tlapchn+(ptau5[jvar][nlevs-k-2][jloc]-ptau5[jvar][nlevs-k-1][jloc])*
                    (tvp[nlevs-k][jloc]-tvp[nlevs-k-2][jloc]);
        }
        tlapchn_all[jvar][jloc] = tlapchn;
      }
    }
  }

  MeansSlot & slot = meansByObsSpace_[&odb];
  if (!slot.ready) {
    // Freeze means after this first fill (file-like lifetime across compute calls).
    slot.means.resize(nvars);
    if (useOnlineMeans_) {
      // First compute on this ObsSpace (typically trajectory / setTrajectory or first H(x)).
      std::unique_ptr<ioda::Accumulator<std::vector<size_t>>> countAcc =
          odb.distribution()->createAccumulator<size_t>(nvars);
      std::unique_ptr<ioda::Accumulator<std::vector<double>>> sumAcc =
          odb.distribution()->createAccumulator<double>(nvars);
      for (std::size_t jloc = 0; jloc < nlocs; ++jloc) {
        for (std::size_t jvar = 0; jvar < nvars; ++jvar) {
          const float tlapchn = tlapchn_all[jvar][jloc];
          if (!isBad(tlapchn)) {
            sumAcc->addTerm(jloc, jvar, static_cast<double>(tlapchn));
            countAcc->addTerm(jloc, jvar, 1);
          }
        }
      }
      const std::vector<size_t> counts = countAcc->computeResult();
      const std::vector<double> sums = sumAcc->computeResult();
      for (std::size_t jvar = 0; jvar < nvars; ++jvar) {
        if (counts[jvar] == 0) {
          oops::Log::error() << "No valid lapse rates to form online mean for channel: "
                             << vars_.channels()[jvar] << std::endl;
          ABORT("No valid lapse rates to form online mean");
        }
        slot.means[jvar] = static_cast<float>(sums[jvar] / static_cast<double>(counts[jvar]));
      }
    } else {
      // sort out the tlapmean based on vars
      for (std::size_t jvar = 0; jvar < nvars; ++jvar) {
        auto it = tlapmean_.find(vars_.channels()[jvar]);
        if (it != tlapmean_.end()) {
          slot.means[jvar] = it->second;
        } else {
          oops::Log::error() << "Could not locate tlapsemean for channel: " <<
                                vars_.channels()[jvar] << std::endl;
          ABORT("Could not locate tlapsemean value");
        }
      }
    }
    slot.ready = true;
  }
  ASSERT(slot.means.size() == nvars);

  for (std::size_t jloc = 0; jloc < nlocs; ++jloc) {
    for (std::size_t jvar = 0; jvar < nvars; ++jvar) {
      const float tlapchn = tlapchn_all[jvar][jloc];
      if (isBad(tlapchn)) {
        out[jloc*nvars+jvar] = dmiss;
      } else {
        out[jloc*nvars+jvar] = pow((tlapchn - slot.means[jvar]), order_);
      }
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
