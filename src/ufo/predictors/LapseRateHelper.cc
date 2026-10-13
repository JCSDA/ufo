/*
 * (C) Copyright 2020 UCAR
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

#include "ufo/predictors/LapseRateHelper.h"

#include "ioda/distribution/Accumulator.h"
#include "ioda/ObsSpace.h"

#include "oops/base/ObsVariables.h"
#include "oops/base/Variables.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

namespace {

bool isBad(float v) {
  const float fmiss = util::missingValue<float>();
  return v == fmiss || std::isnan(v);
}

float missingFloat() {
  return util::missingValue<float>();
}

// Shared per-channel mean lapse rates for all LapseRateHelper instances on one ObsSpace.
// One static map per simulation is keyed to unique ObsSpace objects. This is typical
// for most real applications, and ensures numerical consistency between NL and TL
// in oops TL/AD tests, where some drivers hold two ObsBias objects for the same ObsSpace.
struct ObsSpaceMeanCache {
  bool ready = false;
  std::vector<float> means;
};
std::map<const ioda::ObsSpace *, ObsSpaceMeanCache> cacheByObsSpace_;

}  // namespace

// -----------------------------------------------------------------------------

oops::Variables LapseRateHelper::requiredGeovars() {
  return oops::Variables({"air_temperature",
                          "air_pressure",
                          "average_surface_temperature_within_field_of_view"});
}

oops::ObsVariables LapseRateHelper::requiredHdiags() const {
  return oops::ObsVariables({"transmittances_of_atmosphere_layer"}, channels_);
}

// -----------------------------------------------------------------------------

LapseRateHelper::LapseRateHelper(const LapseRateMeanParameters & parameters,
                                 const oops::ObsVariables & vars)
  : onlineMean_(parameters.onlineMean.value()),
    saveFile_(parameters.saveFile.value().value_or("")) {

  channels_.assign(vars.channels().begin(), vars.channels().end());

  const bool hasTlapseFile = parameters.tlapse.value() != boost::none;
  if (hasTlapseFile == onlineMean_) {
    throw eckit::UserError(
      "LapseRate predictors must specify exactly one of 'tlapse' (file path) or "
      "'online mean: true'",
      Here());
  }

  if (!onlineMean_) {
    // Offline file reading
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
        fileMeans_[nuchan] = tlapse;
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

void LapseRateHelper::fillTlapchn(const ioda::ObsSpace & odb,
                                  const GeoVaLs & geovals,
                                  const ObsDiagnostics & ydiags,
                                  std::vector<std::vector<float>> & tlapchn_by_chan) const {
  const std::size_t nchannels = channels_.size();
  const std::size_t nlocs = odb.nlocs();

  // common vectors storage
  std::vector<float> pred(nlocs, 0.0f);
  // retrieve the average surface temperature
  std::vector<float> tsavg5(nlocs, 0.0f);
  geovals.get(tsavg5, oops::Variable{"average_surface_temperature_within_field_of_view"});

  // Retrieve the transmittances_of_atmosphere_layer from Hdiag
  std::vector<std::vector<std::vector<float>>> ptau5;
  std::vector<std::vector<float>> tmpvar;

  std::string hdiags;
  for (std::size_t ich = 0; ich < nchannels; ++ich) {
    hdiags = "transmittances_of_atmosphere_layer_" + std::to_string(channels_[ich]);
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

  tlapchn_by_chan.assign(nchannels, std::vector<float>(nlocs, missingFloat()));
  if (nchannels == 0 || ptau5.empty() || ptau5[0].size() < nlevs || tvp.size() < nlevs ||
      nlevs < 2) {
    oops::Log::error() << "No valid transmittance profiles or air temperature profiles "
                       << "to compute lapse rates" << std::endl;
    ABORT("No valid transmittance profiles or air temperature profiles to compute lapse rates");
  }

  float tlapchn;
  for (std::size_t ich = 0; ich < nchannels; ++ich) {
    for (std::size_t jloc = 0; jloc < nlocs; ++jloc) {
      // Check for missing or NaN values.
      bool missingInput = isBad(tsavg5[jloc]);
      for (std::size_t lev = 0; lev < nlevs && !missingInput; ++lev) {
        missingInput = isBad(tvp[lev][jloc]) ||
                       isBad(ptau5[ich][lev][jloc]);
      }
      if (missingInput) {
        continue;
      }

      tlapchn = (ptau5[ich][nlevs-2][jloc]-ptau5[ich][nlevs-1][jloc])*
                (tsavg5[jloc]-tvp[nlevs-2][jloc]);
      for (std::size_t k = 1; k < nlevs-1; ++k) {
        tlapchn = tlapchn+(ptau5[ich][nlevs-k-2][jloc]-ptau5[ich][nlevs-k-1][jloc])*
                  (tvp[nlevs-k][jloc]-tvp[nlevs-k-2][jloc]);
      }
      tlapchn_by_chan[ich][jloc] = tlapchn;
    }
  }
}

// -----------------------------------------------------------------------------

void LapseRateHelper::fillTlapchnAndMeans(
    const ioda::ObsSpace & odb,
    const GeoVaLs & geovals,
    const ObsDiagnostics & ydiags,
    std::vector<std::vector<float>> & tlapchn_by_chan,
    std::vector<float> & means_by_chan) const {

  fillTlapchn(odb, geovals, ydiags, tlapchn_by_chan);

  ObsSpaceMeanCache & cache = cacheByObsSpace_[&odb];
  const std::size_t nchannels = channels_.size();

  if (odb.nchans() != nchannels) {
    oops::Log::error() << "Channel mismatch between ObsSpace / LapseRateHelper vars" << std::endl;
    ABORT("Channel mismatch between ObsSpace / LapseRateHelper vars");
  }

  if (!cache.ready) {
    // Freeze means after this first fill (file-like lifetime across fillTlapchnAndMeans calls).
    cache.means.resize(nchannels);
    if (onlineMean_) {
      // First compute on this ObsSpace (typically trajectory / setTrajectory or first H(x)).
      const std::size_t nlocs = odb.nlocs();
      std::unique_ptr<ioda::Accumulator<std::vector<size_t>>> countAcc =
          odb.distribution()->createAccumulator<size_t>(nchannels);
      std::unique_ptr<ioda::Accumulator<std::vector<double>>> sumAcc =
          odb.distribution()->createAccumulator<double>(nchannels);
      for (std::size_t jloc = 0; jloc < nlocs; ++jloc) {
        for (std::size_t ich = 0; ich < nchannels; ++ich) {
          const float v = tlapchn_by_chan[ich][jloc];
          if (!isBad(v)) {
            sumAcc->addTerm(jloc, ich, static_cast<double>(v));
            countAcc->addTerm(jloc, ich, 1);
          }
        }
      }
      const std::vector<size_t> counts = countAcc->computeResult();
      const std::vector<double> sums = sumAcc->computeResult();
      for (std::size_t ich = 0; ich < nchannels; ++ich) {
        if (counts[ich] == 0) {
          cache.means[ich] = missingFloat();
        } else {
          cache.means[ich] =
              static_cast<float>(sums[ich] / static_cast<double>(counts[ich]));
        }
      }
    } else {
      // sort out the tlapmean based on vars
      for (std::size_t ich = 0; ich < nchannels; ++ich) {
        auto it = fileMeans_.find(channels_[ich]);
        if (it != fileMeans_.end()) {
          cache.means[ich] = it->second;
        } else {
          oops::Log::error() << "Could not locate tlapsemean for channel: "
                             << channels_[ich] << std::endl;
          ABORT("Could not locate tlapsemean value");
        }
      }
    }
    cache.ready = true;

    // optionally save means to text file
    if (!saveFile_.empty()) {
      std::ofstream outfile(saveFile_);
      for (std::size_t ich = 0; ich < nchannels; ++ich) {
        outfile << odb.obsname() << " " << channels_[ich] << " " << cache.means[ich] << std::endl;
      }
      outfile.close();
    }
  } else {
    if (!saveFile_.empty()) {
      oops::Log::warning() <<
        "`save means` is only supported on the first lapseRate or cloudWaterContent predictor,"
        " and ignored on all others" << std::endl;
    }
  }

  ASSERT(cache.means.size() == nchannels);
  means_by_chan = cache.means;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
