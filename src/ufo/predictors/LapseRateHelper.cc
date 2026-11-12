/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"

#include "ufo/predictors/LapseRateHelper.h"

#include "ioda/distribution/Accumulator.h"
#include "ioda/distribution/Distribution.h"
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
  bool calculatedMean = false;
  std::string tlapseFile;
  std::string savedFiles;
  std::vector<float> means;
};
std::map<const ioda::ObsSpace *, ObsSpaceMeanCache> cacheByObsSpace_;

// Writes per-channel means to a text file (readable back through `tlapse`) on MPI rank 0,
// overwriting any existing file.
void writeMeans(const ioda::ObsSpace & odb, const std::string & path,
                const std::vector<int> & channels, const std::vector<float> & means) {
  if (odb.comm().rank() != 0) return;
  std::string sensorName = odb.obsname();
  std::replace(sensorName.begin(), sensorName.end(), ' ', '_');
  std::ofstream outfile(path);
  if (!outfile.is_open()) {
    oops::Log::error() << "LapseRateHelper: Unable to open file : " << path << std::endl;
    ABORT("LapseRateHelper: Unable to open save means file ");
  }
  for (std::size_t ich = 0; ich < channels.size(); ++ich) {
    outfile << sensorName << " " << channels[ich] << " "
            << std::setprecision(std::numeric_limits<float>::max_digits10)
            << means[ich] << "\n";
  }
}

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
  : calculatedMean_(parameters.calculatedMean.value()),
    tlapseFile_(parameters.tlapseFile.value().value_or("")),
    saveFile_(parameters.saveFile.value().value_or("")) {
  channels_.assign(vars.channels().begin(), vars.channels().end());

  const bool hasTlapseFile = !tlapseFile_.empty();
  if (hasTlapseFile == calculatedMean_) {
    throw eckit::UserError(
      "LapseRate predictors must specify exactly one of 'tlapse' (file path) or "
      "'calculated mean: true'",
      Here());
  }

  if (hasTlapseFile) {
    readTlapseFile();
  }
}

// -----------------------------------------------------------------------------

void LapseRateHelper::readTlapseFile() {
  // Each non-blank line is "<sensor> <channel> <mean>". The sensor name is not used and may
  // contain spaces, so the channel and mean are taken from the last two fields.
  std::ifstream infile(tlapseFile_);
  if (!infile.is_open()) {
    oops::Log::error() << "Unable to open file : " << tlapseFile_ << std::endl;
    ABORT("Unable to open tlapse file ");
  }

  std::string line;
  std::size_t lineNumber = 0;
  while (std::getline(infile, line)) {
    ++lineNumber;
    std::istringstream iss(line);
    std::vector<std::string> fields;
    std::string field;
    while (iss >> field) fields.push_back(field);
    if (fields.empty()) continue;

    bool goodLine = fields.size() >= 3;
    int nuchan = 0;
    float tlapse = 0.0f;
    if (goodLine) {
      std::istringstream chanStream(fields[fields.size() - 2]);
      std::istringstream meanStream(fields[fields.size() - 1]);
      goodLine = (chanStream >> nuchan) && chanStream.eof() &&
               (meanStream >> tlapse) && meanStream.eof();
    }
    if (!goodLine) {
      oops::Log::error() << "LapseRateHelper: malformed line " << lineNumber << " in tlapse file "
                         << tlapseFile_ << ": '" << line << "'" << std::endl;
      ABORT("LapseRateHelper: malformed tlapse file");
    }
    fileMeans_[nuchan] = tlapse;
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

  tlapchn_by_chan.assign(nchannels, std::vector<float>(nlocs, missingFloat()));
  const bool transmittanceTooShort =
      std::any_of(ptau5.begin(), ptau5.end(),
                  [nlevs](const std::vector<std::vector<float>> & p) { return p.size() < nlevs; });
  if (nchannels == 0 || transmittanceTooShort || nlevs < 2) {
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
        continue;  // tlapchn_by_chan[ich][jloc] stays missing
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

  if (!cache.ready) {
    // Freeze means after this first fill (file-like lifetime across fillTlapchnAndMeans calls).
    cache.calculatedMean = calculatedMean_;
    cache.tlapseFile = tlapseFile_;
    cache.means.assign(nchannels, missingFloat());
    if (calculatedMean_) {
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
        if (counts[ich] > 0) {
          cache.means[ich] = static_cast<float>(sums[ich] / static_cast<double>(counts[ich]));
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
          ABORT("LapseRateHelper: Could not locate tlapsemean value");
        }
      }
    }
    cache.ready = true;
  } else {
    if (cache.calculatedMean != calculatedMean_) {
      oops::Log::error() << "LapseRateHelper: inconsistent mean source for ObsSpace '"
                         << odb.obsname() << "': cache uses "
                         << (cache.calculatedMean ? "'calculated mean'" : "'tlapse' file")
                         << " but this predictor uses "
                         << (calculatedMean_ ? "'calculated mean'" : "'tlapse' file")
                         << std::endl;
      ABORT("LapseRateHelper: inconsistent calculated mean / tlapse configuration");
    }
    if (!calculatedMean_ && cache.tlapseFile != tlapseFile_) {
      oops::Log::error() << "LapseRateHelper: inconsistent tlapse file for ObsSpace '"
                         << odb.obsname() << "': cache uses '" << cache.tlapseFile
                         << "' but this predictor uses '" << tlapseFile_ << "'"
                         << std::endl;
      ABORT("LapseRateHelper: inconsistent tlapse file configuration");
    }
  }

  ASSERT(cache.means.size() == nchannels);

  // optionally save the frozen means; any predictor on this ObsSpace may request a path,
  // and the file is written once
  if (!saveFile_.empty()) {
    if (cache.savedFiles.empty()) {
      writeMeans(odb, saveFile_, channels_, cache.means);
      cache.savedFiles = saveFile_;
    } else if (cache.savedFiles != saveFile_) {
      oops::Log::error() << "LapseRateHelper: inconsistent save means file for ObsSpace '"
                         << odb.obsname() << "': means already saved to '" << cache.savedFiles
                         << "' but this predictor requests '" << saveFile_ << "'"
                         << std::endl;
      ABORT("LapseRateHelper: inconsistent save means file configuration");
    }
  }
  means_by_chan = cache.means;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
