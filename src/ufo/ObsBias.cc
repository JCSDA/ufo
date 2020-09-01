/*
 * (C) Copyright 2017-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/ObsBias.h"

#include <Eigen/Core>
#include <fstream>
#include <iomanip>
#include <set>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"

#include "ufo/ObsBiasIncrement.h"
#include "ufo/ObsDiagnostics.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsBias::ObsBias(ioda::ObsSpace & odb, const eckit::Configuration & conf)
  : predbases_(0), jobs_(0), odb_(odb), conf_(conf) {
  oops::Log::trace() << "ObsBias::create starting." << std::endl;

  // Get the jobs(channels)
  if (conf_.has("obs bias.jobs")) {
    const std::set<int> jobs = oops::parseIntSet(conf_.getString("obs bias.jobs"));
    jobs_.assign(jobs.begin(), jobs.end());
  }

  // Predictor factory
  if (conf_.has("obs bias.predictors")) {
    std::vector<eckit::LocalConfiguration> confs;
    conf_.get("obs bias.predictors", confs);
    for (std::size_t j = 0; j < confs.size(); ++j) {
      std::shared_ptr<PredictorBase> pred(PredictorFactory::create(confs[j], jobs_));
      predbases_.push_back(pred);
      prednames_.push_back(pred->name());
      geovars_ += pred->requiredGeovars();
      hdiags_ += pred->requiredHdiagnostics();

      // Reserve the space for ObsBiasTerm for predictor
      if (jobs_.size() > 0) {
        for (auto & job : jobs_) {
          hdiags_ += oops::Variables({prednames_.back() + "_" + std::to_string(job)});
        }
      } else {
        hdiags_ += oops::Variables({prednames_.back()});
      }
    }
  }

  if (prednames_.size() * jobs_.size() > 0) {
    // Initialize the biascoeffs to ZERO
    biascoeffs_.resize(prednames_.size() * jobs_.size(), 0.0);

    // Read or initialize bias coefficients
    this->read(conf);
  }

  oops::Log::trace() << "ObsBias::create done." << std::endl;
}

// -----------------------------------------------------------------------------

ObsBias::ObsBias(const ObsBias & other, const bool copy)
  : odb_(other.odb_), conf_(other.conf_), predbases_(other.predbases_),
    prednames_(other.prednames_), jobs_(other.jobs_),
    geovars_(other.geovars_), hdiags_(other.hdiags_) {
  oops::Log::trace() << "ObsBias::copy ctor starting." << std::endl;

  // Initialize the biascoeffs
  biascoeffs_.resize(prednames_.size() * jobs_.size(), 0.0);

  // Copy the bias coeff data
  if (copy && biascoeffs_.size() > 0) *this = other;

  oops::Log::trace() << "ObsBias::copy ctor done." << std::endl;
}

// -----------------------------------------------------------------------------

ObsBias & ObsBias::operator+=(const ObsBiasIncrement & dx) {
  for (std::size_t jj = 0; jj < biascoeffs_.size(); ++jj)
      biascoeffs_[jj] += dx[jj];
  return *this;
}

// -----------------------------------------------------------------------------

ObsBias & ObsBias::operator=(const ObsBias & rhs) {
  if (rhs.size() > 0 && this->size() == rhs.size()) {
    conf_       = rhs.conf_;
    biascoeffs_ = rhs.biascoeffs_;
    predbases_  = rhs.predbases_;
    prednames_  = rhs.prednames_;
    jobs_       = rhs.jobs_;
    geovars_    = rhs.geovars_;
    hdiags_     = rhs.hdiags_;
  }
  return *this;
}

// -----------------------------------------------------------------------------

void ObsBias::read(const eckit::Configuration & conf) {
  oops::Log::trace() << "ObsBias::read and initialize from file, starting "<< std::endl;

  // Bias coefficients input file name
  std::string input_filename;
  if (conf.has("obs bias.abias_in")) {
    input_filename = conf.getString("obs bias.abias_in");
    oops::Log::debug() << "ObsBias::initialize coefficients from file: "
                       << input_filename <<  std::endl;

    // Default predictor names from GSI
    // temporary solution, we should have a self-explanatory obsbias file
    const std::vector<std::string> gsi_predictors = {"constant",
                                                     "zenith_angle",
                                                     "cloud_liquid_water",
                                                     "lapse_rate_order_2",
                                                     "lapse_rate",
                                                     "cosine_of_latitude_times_orbit_node",
                                                     "sine_of_latitude",
                                                     "emissivity",
                                                     "scan_angle_order_4",
                                                     "scan_angle_order_3",
                                                     "scan_angle_order_2",
                                                     "scan_angle"
                                                     };
    std::ifstream infile(input_filename);

    // get the sensor id
    std::string sensor = conf.getString("obs bias.sensor");
    oops::Log::debug() << "ObsBias::initialize coefficients for sensor: "
                       << sensor <<  std::endl;

    std::size_t ich;     //  sequential number
    std::string nusis;   //  sensor/instrument/satellite
    std::size_t nuchan;  //  channel number
    float tlap, tsum;
    std::size_t ntlapupdate;

    if (infile.is_open())
    {
      oops::Log::debug() << "ObsBias:: prior file is opened" << std::endl;
      float par;
      while (!infile.eof())
      {
        infile >> ich;
        infile >> nusis;
        infile >> nuchan;
        infile >> tlap;
        infile >> tsum;
        infile >> ntlapupdate;
        if (nusis == sensor) {
          auto ijob = std::find(jobs_.begin(), jobs_.end(), nuchan);
          if (ijob != jobs_.end()) {
            int j = std::distance(jobs_.begin(), ijob);

            for (auto & item : gsi_predictors) {
              infile >> par;
              auto ipred = std::find(prednames_.begin(), prednames_.end(), item);
              if (ipred != prednames_.end()) {
                int p = std::distance(prednames_.begin(), ipred);
                biascoeffs_.at(j*prednames_.size() + p) = static_cast<double>(par);
              }
            }
          }
        } else {
          for (auto item : gsi_predictors) {
            infile >> par;
          }
        }
      }
      infile.close();
      oops::Log::debug() << "ObsBias:: read prior from " << input_filename << std::endl;
    } else {
      oops::Log::error() << "Unable to open file : " << input_filename << std::endl;
      ABORT("Unable to open bias prior file ");
    }
  } else {
    oops::Log::warning() << "ObsBias::prior file is NOT available, starting from ZERO"
                         << std::endl;
  }

  oops::Log::trace() << "ObsBias::read and initilization done " << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBias::write(const eckit::Configuration & conf) const {
  oops::Log::trace() << "ObsBias::write to file not implemented" << std::endl;

  // Bias coefficients output file name
  std::string output_filename;
  if (conf.has("ObsBias.abias_out")) {
    output_filename = conf.getString("ObsBias.abias_out");
  }

  oops::Log::trace() << "ObsBias::write to file done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBias::computeObsBias(ioda::ObsVector & ybias,
                             ObsDiagnostics & ydiags,
                             const ioda::ObsDataVector<double> & predData) const {
  oops::Log::trace() << "ObsBias::computeObsBias starting" << std::endl;

  if (this->size() > 0) {
    const std::size_t nlocs  = ybias.nlocs();
    const std::size_t npreds = prednames_.size();
    const std::size_t njobs  = jobs_.size();

    ASSERT(biascoeffs_.size() == npreds*njobs);
    ASSERT(predData.nlocs() == nlocs);
    ASSERT(predData.nvars() == njobs*npreds);
    ASSERT(ybias.nvars() == njobs);

    /* predData memory layout (njobs*npreds X nlocs)
     *     Loc     0      1      2       3
     *           --------------------------
     * ch1 pred1 | 0      1      2       3
     *     pred2 | 4      5      6       7
     *     pred3 | 8      9     10      11 
     * ch2 pred1 |12     13     14      15
     *     pred2 |16     17     18      19
     *     ....  |
     */

    ybias.zero();

    /* ybias memory layout (nlocs X njobs)
     *     ch1    ch2    ch3     ch4
     * Loc --------------------------
     *  0 | 0      1      2       3
     *  1 | 4      5      6       7
     *  2 | 8      9     10      11 
     *  3 |12     13     14      15
     *  4 |16     17     18      19
     * ...|
     */

    /* map bias coeff to eigen matrix npreds X njobs (read only)
     * bias coeff memory layout (npreds X njobs)
     *        ch1    ch2    ch3     ch4
     *       --------------------------
     * pred1 | 0      1      2       3
     * pred2 | 4      5      6       7
     * pred3 | 8      9     10      11 
     * ....  |
     */
    Eigen::Map<const Eigen::MatrixXd> coeffs(biascoeffs_.data(), npreds, njobs);

    std::vector<double> biasTerm(nlocs, 0.0);
    //  ( nlocs X 1 ) =  ( nlocs X npreds ) * (  npreds X 1 )
    for (std::size_t jch = 0; jch < njobs; ++jch) {
      for (std::size_t jp = 0; jp < npreds; ++jp) {
        const std::string varname = predbases_[jp]->name() + "_" + std::to_string(jobs_[jch]);
        for (std::size_t jl = 0; jl < nlocs; ++jl) {
          biasTerm[jl] = predData[varname][jl] * coeffs(jp, jch);
          ybias[jl*njobs+jch] += biasTerm[jl];
        }
        // Save ObsBiasTerms (bias_coeff x predictor) for QC
        if (ydiags.has(varname)) {
          ydiags.save(biasTerm, varname, 1);
        } else {
          oops::Log::error() << varname << " is not reserved in ydiags !" << std::endl;
          ABORT("ObsBiasTerm variable is not reserved in ydiags");
        }
      }
    }
  }

  oops::Log::trace() << "ObsBias::computeObsBias done." << std::endl;
}

// -----------------------------------------------------------------------------
ioda::ObsDataVector<double> ObsBias::computePredictors(const GeoVaLs & geovals,
                                                       const ObsDiagnostics & ydiags) const {
  const std::size_t nlocs  = odb_.nlocs();
  const std::size_t npreds = predbases_.size();
  const std::size_t njobs  = jobs_.size();

  const oops::Variables vars(prednames_, jobs_);
  ioda::ObsDataVector<double> predData(odb_, vars);

  for (std::size_t r = 0; r < npreds; ++r) {
    predbases_[r]->compute(odb_, geovals, ydiags, predData);
  }

  oops::Log::trace() << "ObsBias::computePredictors done." << std::endl;
  return predData;
}

// -----------------------------------------------------------------------------

double ObsBias::norm() const {
  oops::Log::trace() << "ObsBias::norm starting." << std::endl;
  double zz = 0.0;
  for (std::size_t jj = 0; jj < biascoeffs_.size(); ++jj) {
    zz += biascoeffs_[jj] * biascoeffs_[jj];
  }
  if (biascoeffs_.size() > 0) zz = std::sqrt(zz/biascoeffs_.size());
  oops::Log::trace() << "ObsBias::norm done." << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

void ObsBias::print(std::ostream & os) const {
  if (this->size() > 0) {
    // map bias coeffs to eigen matrix (writable)
    Eigen::Map<const Eigen::MatrixXd>
      coeffs(biascoeffs_.data(), prednames_.size(), jobs_.size());

    os << "ObsBias::print " << std::endl;
    os << "---------------------------------------------------------------" << std::endl;
    auto njobs = jobs_.size();
    for (std::size_t p = 0; p < prednames_.size(); ++p) {
      os << std::fixed << std::setw(20) << prednames_[p]
         << ":  Min= " << std::setw(15) << std::setprecision(8)
         << coeffs.row(p).minCoeff()
         << ",  Max= " << std::setw(15) << std::setprecision(8)
         << coeffs.row(p).maxCoeff()
         << ",  Norm= " << std::setw(15) << std::setprecision(8)
         << coeffs.row(p).norm()
         << std::endl;
    }
    os << "---------------------------------------------------------------" << std::endl;
    os << "ObsBias::print done" << std::endl;
  }
}

// -----------------------------------------------------------------------------

}  // namespace ufo
