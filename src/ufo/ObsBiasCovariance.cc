/*
 * (C) Copyright 2018-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <memory>
#include <random>
#include <set>

#include "ufo/ObsBiasCovariance.h"

#include "ioda/distribution/Accumulator.h"
#include "ioda/Engines/EngineUtils.h"
#include "ioda/Engines/HH.h"
#include "ioda/Layout.h"
#include "ioda/ObsGroup.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"

#include "oops/util/IntSetParser.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "oops/util/Random.h"

#include "ufo/ObsBias.h"
#include "ufo/ObsBiasIncrement.h"
#include "ufo/ObsBiasPreconditioner.h"
#include "ufo/predictors/PredictorBase.h"
#include "ufo/utils/IodaGroupIndices.h"

namespace ufo {

// -----------------------------------------------------------------------------

ObsBiasCovariance::ObsBiasCovariance(ioda::ObsSpace & odb, const eckit::Configuration & config)
  : odb_(odb), ht_rinv_h_(0), preconditioner_(0), obs_num_(0),
    minimal_required_obs_number_(0), analysis_variances_(0), variances_(),
    prednames_(0), vars_(odb.assimvariables()), rank_(odb.distribution()->rank()),
    commTime_(odb.commTime())
{
  oops::Log::trace() << "ObsBiasCovariance::Constructor starting" << std::endl;

  Parameters_ params;
  params.validateAndDeserialize(config);
  byRecord_ = params.BiasCorrectionByRecord;

  // Predictor factory
  for (const PredictorParametersWrapper &wrapper :
       params.variationalBC.value().predictors.value()) {
    std::shared_ptr<PredictorBase> pred(PredictorFactory::create(wrapper.predictorParameters,
                                                                 vars_));
    prednames_.push_back(pred->name());
  }
  nrecs_ = (byRecord_ && odb.obs_group_vars().size() > 0) ? odb.nrecs() : 1;
  if (byRecord_ && odb.obs_group_vars().size() == 0) {
    throw eckit::BadParameter("ObsBiasParameters: BiasCorrectionByRecord is turned on, "
                              "but the observations are not grouped into records.");
  }
  ASSERT(nrecs_ > 0);

  oops::ObsVariables varsNoBC = params.variablesNoBC;
  varsNoBC.intersection(vars_);  // Safeguard to make sure that varsNoBC is a subset of vars_
  for (size_t ii = 0; ii < varsNoBC.size(); ++ii) {
    size_t index = vars_.find(varsNoBC[ii]);
    varIndexNoBC_.push_back(index);
  }

  // save record IDs for matching
  if (byRecord_) {
    recIds_.resize(nrecs_);
    // get all ids and obs types (for the hack to be removed)
    std::vector<std::string> allids;
    odb.get_db("MetaData", "stationIdentification", allids);
    // save station ids for all records
    size_t jrec = 0;
    for (auto irec = odb.recidx_begin(); irec != odb.recidx_end(); ++irec, ++jrec) {
      // all the identifiers will be the same for the same record, use the first one
      const size_t iloc = odb.recidx_vector(irec)[0];
      // remove trailing whitespaces (should really be done in files)
      const size_t strEnd = allids[iloc].find_last_not_of(" \t");
      recIds_[jrec] = allids[iloc].substr(0, strEnd+1);
    }
  }

  if (vars_.size() * prednames_.size() > 0) {
    if (params.covariance.value() == boost::none)
      throw eckit::UserError("obs bias.covariance section missing from the YAML file");
    const ObsBiasCovarianceParameters &biasCovParams = *params.covariance.value();

    // Get the minimal required filtered obs number
    minimal_required_obs_number_ = biasCovParams.minimalRequiredObsNumber;

    // Override the variance range if provided
    {
      const std::vector<double> &range = biasCovParams.varianceRange.value();
      ASSERT(range.size() == 2);
      smallest_variance_ = range[0];
      largest_variance_ = range[1];
    }

    // Override the preconditioning step size if provided
    step_size_ = biasCovParams.stepSize;

    // Override the largest analysis variance if provided
    largest_analysis_variance_ = biasCovParams.largestAnalysisVariance;

    // Set variance for new record ("new record analysis variance" if specified
    // in yaml, or largest_analysis_variance if not))
    default_new_variance_ = (biasCovParams.defaultNewVariance.value() != boost::none) ?
            *biasCovParams.defaultNewVariance.value() : largest_analysis_variance_;

    // Initialize the variances to upper limit
    variances_ = Eigen::VectorXd::Constant(nrecs_ * vars_.size() * prednames_.size(),
                                           largest_variance_);

    // Initialize the hessian contribution to zero
    ht_rinv_h_.resize(nrecs_ * vars_.size() * prednames_.size());
    std::fill(ht_rinv_h_.begin(), ht_rinv_h_.end(), 0.0);

    // Initialize the preconditioner to default step size
    preconditioner_.resize(nrecs_ * vars_.size() * prednames_.size());
    std::fill(preconditioner_.begin(), preconditioner_.end(), step_size_);

    // Initialize obs_num_ to ZERO
    obs_num_.resize(nrecs_ * vars_.size());
    std::fill(obs_num_.begin(), obs_num_.end(), 0);

    // Initialize analysis error variances to the upper limit
    analysis_variances_.resize(nrecs_ * vars_.size() * prednames_.size());
    std::fill(analysis_variances_.begin(), analysis_variances_.end(), largest_analysis_variance_);

    // Initializes from given prior
    if (biasCovParams.prior.value() != boost::none) {
      const ObsBiasCovariancePriorParameters &priorParams = *biasCovParams.prior.value();

      // Get default inflation ratio
      const double inflation_ratio = priorParams.inflation.value().ratio;

      // Check the large inflation ratio when obs number < minimal_required_obs_number
      const double large_inflation_ratio = priorParams.inflation.value().ratioForSmallDataset;

      // read in Variances prior (analysis_variances_) and number of obs. (obs_num_)
      // from previous cycle
      this->read(priorParams.toConfiguration());

      // set variances for bias predictor coeff. based on diagonal info
      // of previous analysis error variance
      std::size_t ii;
      for (std::size_t j = 0; j < nrecs_ * vars_.size(); ++j) {
        const double inflation = (obs_num_[j] <= minimal_required_obs_number_) ?
                                 large_inflation_ratio : inflation_ratio;
        for (std::size_t p = 0; p < prednames_.size(); ++p) {
          ii = j * prednames_.size() + p;
          if (inflation > inflation_ratio) {
             analysis_variances_[ii] = inflation * analysis_variances_[ii] + smallest_variance_;
             variances_[ii] = analysis_variances_[ii];
          } else {
             variances_[ii] = inflation_ratio * analysis_variances_[ii] + smallest_variance_;
          }
          if (variances_[ii] > largest_variance_) variances_[ii] = largest_variance_;
          if (analysis_variances_[ii] > largest_analysis_variance_)
            analysis_variances_[ii] = largest_analysis_variance_;
        }
      }
    }
  }

  oops::Log::trace() << "ObsBiasCovariance::Constructor is done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasCovariance::read(const eckit::Configuration & config) {
  oops::Log::trace() << "ObsBiasCovariance::read from file " << std::endl;

  ObsBiasCovariancePriorParameters params;
  params.validateAndDeserialize(config);

  if (params.inputFile.value() != boost::none) {
    // Open an hdf5 file, read only
    ioda::Engines::BackendNames  backendName = ioda::Engines::BackendNames::Hdf5File;
    ioda::Engines::BackendCreationParameters backendParams;
    backendParams.fileName = *params.inputFile.value();
    backendParams.action   = ioda::Engines::BackendFileActions::Open;
    backendParams.openMode = ioda::Engines::BackendOpenModes::Read_Only;

    // Create the backend and attach it to an ObsGroup
    // Use the None DataLyoutPolicy for now to accommodate the current file format
    ioda::Group backend = constructBackend(backendName, backendParams);
    ioda::ObsGroup obsgroup = ioda::ObsGroup(backend,
                   ioda::detail::DataLayoutPolicy::generate(
                         ioda::detail::DataLayoutPolicy::Policies::None));

    // setup variables
    std::vector<Eigen::ArrayXXf> allbcerrors;
    std::vector<std::string> predictors;

    // loop through list of coefficients, read them, and store in vector
    for (size_t jpred = 0; jpred < prednames_.size(); ++jpred) {
      ioda::Variable bcerrvar = obsgroup.vars["BiasCoefficientErrors/"+prednames_[jpred]];
      Eigen::ArrayXXf bcerrs;
      bcerrvar.readWithEigenRegular(bcerrs);
      allbcerrors.push_back(bcerrs);
      predictors.push_back(prednames_[jpred]);
    }

    // Read nobs into Eigen array
    ioda::Variable nobsvar = obsgroup.vars["numberObservationsUsed"];
    Eigen::ArrayXXf nobsassim;
    nobsvar.readWithEigenRegular(nobsassim);

    // Read all record names into the Eigen array
    std::vector<std::string> allrecords;
    if (obsgroup.vars.exists("stationIdentification")) {
      ioda::Variable recvar = obsgroup.vars.open("stationIdentification");
      recvar.read<std::string>(allrecords);
    }

    // Find indices of variables/channels that we need in the data read from the file
    // Don't throw an exception if the variable is not in the file if it does not need to be
    // bias-corrected.
    bool throwexception = (varIndexNoBC_.size() == 0) ? true : false;
    const std::vector<int> var_idx = getRequiredVarOrChannelIndices(obsgroup, vars_,
                                     throwexception);
    // sanity check
    for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
      if (var_idx[jvar] == -1) {
        ASSERT(std::find(varIndexNoBC_.begin(), varIndexNoBC_.end(), jvar) != varIndexNoBC_.end());
      }
    }
    // Find indices of predictors that we need in the data read from the file
    const std::vector<int> pred_idx = getAllStrIndices(predictors,
                                      prednames_.begin(), prednames_.end());
    // Determine if the records are in the input file, if not, add it to the list
    std::vector<int> rec_idx;
    if (byRecord_) {
      bool throwexception = false;
      rec_idx = getAllStrIndices(allrecords,
                recIds_.begin(), recIds_.end(), throwexception);
    } else {
      rec_idx.push_back(0);
    }

    // Filter predictors and channels that we need
    // FIXME: may be possible by indexing allbcerrors(pred_idx, chan_idx) when Eigen 3.4
    // is available
    const double missing = util::missingValue<double>();
    for (size_t jrec = 0; jrec < nrecs_; ++jrec) {
      for (size_t jvar = 0; jvar < var_idx.size(); ++jvar) {
        const size_t jrecvar = jrec * vars_.size() + jvar;
        if (rec_idx[jrec] == -1) {
          obs_num_[jrecvar] = 0;
        } else {
          obs_num_[jrecvar] = nobsassim(rec_idx[jrec], var_idx[jvar]);
        }
        if (rec_idx[jrec] == -1) {
          // errors are set to default new variance if record not in input file
          for (size_t jpred = 0; jpred < pred_idx.size(); ++jpred) {
            analysis_variances_[index(jrec, jvar, jpred)] = default_new_variance_;
          }
        } else if (var_idx[jvar] == -1) {
          // errors are set to missing values if variable is not bias corrected
          for (size_t jpred = 0; jpred < pred_idx.size(); ++jpred) {
            analysis_variances_[index(jrec, jvar, jpred)] = missing;
          }
        } else {
          for (size_t jpred = 0; jpred < pred_idx.size(); ++jpred) {
            analysis_variances_[index(jrec, jvar, jpred)] =
                 allbcerrors[pred_idx[jpred]](rec_idx[jrec], var_idx[jvar]);
          }
        }
      }
    }
  }
  oops::Log::trace() << "ObsBiasCovariance::read is done " << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasCovariance::write(const eckit::Configuration & config) {
  // only write files out on the task with MPI rank 0
  if (rank_ != 0 || commTime_.rank() != 0) return;

  oops::Log::trace() << "ObsBiasCovariance::write to file " << std::endl;

  Parameters_ params;
  params.validateAndDeserialize(config);

  // only write files out if the obs bias covariance section is present.
  if (params.covariance.value() == boost::none) {
    oops::Log::trace() << "ObsBiasCovariance::write: exiting because the "
                       << "`obs bias.covariance` yaml section is not present."
                       << std::endl;
    return;
  }

  const ObsBiasCovarianceParameters &biasCovParams = *params.covariance.value();

  if (config.has("covariance") && biasCovParams.outputFile.value() != boost::none) {
    // FIXME: only implemented for channels currently
    if (vars_.channels().size() == 0) {
      throw eckit::NotImplemented("ObsBiasCovariance::write not implemented for without channels",
                                  Here());
    }
    // Create a file, overwrite if exists
    const std::string output_filename = *biasCovParams.outputFile.value();
    ioda::Group group = ioda::Engines::HH::createFile(output_filename,
                        ioda::Engines::BackendCreateModes::Truncate_If_Exists);

    ioda::ObsGroup ogrp;

    // put only variable bias predictors into the predictors vector
    std::vector<std::string> predictors(prednames_.begin(), prednames_.end());

    std::vector<int> obs_assimilated(obs_num_.begin(), obs_num_.end());
    Eigen::Map<const Eigen::MatrixXd>
       allbcerrors(analysis_variances_.data(), prednames_.size(), nrecs_ * vars_.size());
    const std::vector<int> channels = vars_.channels();
    // dimensions
    ioda::NewDimensionScales_t dims {
          ioda::NewDimensionScale<int>("Record", 1),
          ioda::NewDimensionScale<int>("Channel", channels.size())
    };
    // new ObsGroup
    ogrp = ioda::ObsGroup::generate(group, dims);
    // and the variables
    ioda::Variable chansVar = ogrp.vars.open("Channel");
    chansVar.write(channels);

    // write number_obs_assimilated
    ioda::Variable nobsVar = ogrp.vars.createWithScales<int>(
                             "numberObservationsUsed", {ogrp.vars["Record"], ogrp.vars["Channel"]});

    nobsVar.write(obs_assimilated);

    // Set up the creation parameters for the bias covariance coefficients variable
    ioda::VariableCreationParameters float_params;
    float_params.chunk = true;               // allow chunking
    float_params.compressWithGZIP();         // compress using gzip
    const float missing_value = util::missingValue<float>();
    float_params.setFillValue<float>(missing_value);

    // Loop over predictors and create variables
    for (size_t jpred = 0; jpred < predictors.size(); ++jpred) {
      // create and write the bias covariance coeffs
      ioda::Variable anvarVar = ogrp.vars.createWithScales<float>(
                               "BiasCoefficientErrors/"+predictors[jpred],
                               {ogrp.vars["Record"], ogrp.vars["Channel"]}, float_params);
      anvarVar.writeWithEigenRegular(allbcerrors(jpred, Eigen::all));
    }
  }
  oops::Log::trace() << "ObsBiasCovariance::write is done " << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasCovariance::linearize(const ObsBias & bias, const eckit::Configuration & innerConf) {
  oops::Log::trace() << "ObsBiasCovariance::linearize starts" << std::endl;
  ASSERT(nrecs_ == bias.nrecs());

  if (vars_.size() * prednames_.size() > 0) {
    const double missing = util::missingValue<double>();
    const int missing_int = util::missingValue<int>();
    const int jouter = innerConf.getInt("iteration");
    // The MPI accumulator only required for not-by-record BC
    std::unique_ptr<ioda::Accumulator<std::vector<size_t>>> obs_num_accumulator =
        odb_.distribution()->createAccumulator<size_t>(obs_num_.size());

    // Retrieve the QC flags of the current outer loop and recalculate the number of
    // effective obs
    const std::string qc_group_name = "EffectiveQC" + std::to_string(jouter);
    const std::vector<std::string> vars = odb_.assimvariables().variables();
    std::vector<int> qc_flags(odb_.nlocs(), 999);
    // get all global records for bias correction by record
    const std::vector<size_t> recnums = odb_.recidx_all_recnums();
    for (std::size_t jvar = 0; jvar < vars.size(); ++jvar) {
      if (odb_.has(qc_group_name, vars[jvar])) {
        odb_.get_db(qc_group_name, vars[jvar], qc_flags);
        for (std::size_t jloc = 0; jloc < qc_flags.size(); ++jloc)
          if (qc_flags[jloc] == 0) {
            size_t jrec = 0;
            if (nrecs_ > 1) {
              std::size_t jrec_global = odb_.recnum()[jloc];
              // find the index for this record on the local task
              jrec = std::find(recnums.begin(), recnums.end(), jrec_global) -
                               recnums.begin();
            }
            const size_t jrecvar = jrec * vars_.size() + jvar;
            if (byRecord_) {
              // For by-record BC each MPI task has information about non-overlapping
              // records, and number of obs can be computed independently
              obs_num_[jrecvar]++;
            } else {
              // For not by-record BC number of obs is cumulative across MPI tasks
              // and requires MPI accumulator
              obs_num_accumulator->addTerm(jloc, jrecvar, 1);
            }
          }
      } else {
        throw eckit::UserError("Unable to find QC flags : " + vars[jvar] + "@" + qc_group_name);
      }
    }

    // For not by-record BC number of obs is cumulative across MPI tasks, sum across
    // the processors
    if (!byRecord_) {
      obs_num_ = obs_num_accumulator->computeResult();
    }
    // Sum across time subwindows
    commTime_.allReduceInPlace(obs_num_.begin(), obs_num_.end(), eckit::mpi::sum());
    // compute the hessian contribution from Jo bias terms channel by channel
    // retrieve the effective error (after QC) for this channel
    const std::string err_group_name = "EffectiveError" + std::to_string(jouter);
    ioda::ObsVector r_inv(odb_, err_group_name);
    // compute \mathrm{R}^{-1}
    std::size_t nvars = r_inv.nvars();
    for (size_t vv = 0; vv < nvars; ++vv) {
      for (size_t ii = 0; ii < r_inv.nlocs(); ++ii) {
        if (r_inv[ii*nvars + vv] != missing) {
          r_inv[ii*nvars + vv] = 1.0f / pow(r_inv[ii*nvars + vv], 2);
        } else {
          r_inv[ii*nvars + vv] = 0.0f;
        }
      }
    }

    // compute \mathrm{H}_\beta^\intercal \mathrm{R}^{-1} \mathrm{H}_\beta
    // -----------------------------------------
    // The MPI accumulator only required for not-by-record BC
    std::unique_ptr<ioda::Accumulator<std::vector<double>>> ht_rinv_h_accumulator =
        odb_.distribution()->createAccumulator<double>(ht_rinv_h_.size());
    for (std::size_t jp = 0; jp < prednames_.size(); ++jp) {
      // retrieve the predictors
      const ioda::ObsVector predx(odb_, prednames_[jp] + "Predictor");
      // for each variable
      ASSERT(r_inv.nlocs() == predx.nlocs());
      // only keep the diagnoal
      for (size_t jvar = 0; jvar < vars_.size(); ++jvar) {
        for (size_t jloc = 0; jloc < predx.nlocs(); ++jloc) {
          size_t jrec = 0;
          if (nrecs_ > 1) {
            std::size_t jrec_global = odb_.recnum()[jloc];
            // find the index for this record on the local task
            jrec = std::find(recnums.begin(), recnums.end(), jrec_global) -
                             recnums.begin();
          }
          if (byRecord_) {
            // For by-record BC each MPI task has information about non-overlapping
            // records, and result can be accumulated independently on each task
            ht_rinv_h_[index(jrec, jvar, jp)] += pow(predx[jloc * vars_.size() + jvar], 2)
                                                  * r_inv[jloc * vars_.size() + jvar];
          } else {
            // For not by-record BC the result is cumulative across MPI tasks
            // and requires MPI accumulator
            ht_rinv_h_accumulator->addTerm(jloc, index(jrec, jvar, jp),
                                           pow(predx[jloc * vars_.size() + jvar], 2)
                                             * r_inv[jloc * vars_.size() + jvar]);
          }
        }
      }
    }

    // For not by-record BC the hessian is cumulative across MPI tasks, sum contributions across
    // the tasks
    if (!byRecord_) {
      ht_rinv_h_ = ht_rinv_h_accumulator->computeResult();
    }
    // Sum across time subwindows
    commTime_.allReduceInPlace(ht_rinv_h_.begin(), ht_rinv_h_.end(), eckit::mpi::sum());
    // Set obs_num_ and ht_rinv_h_ to missing for variables opted out of bias correction
    for (std::size_t jrec = 0; jrec < nrecs_; ++jrec) {
      for (const int jvar : varIndexNoBC_) {
        const std::size_t jrecvar = jrec * vars_.size() + jvar;
        obs_num_[jrecvar] = missing_int;
        for (std::size_t jp = 0; jp < prednames_.size(); ++jp) {
          ht_rinv_h_[index(jrec, jvar, jp)] = missing;
        }
      }
    }

    for (std::size_t jrec = 0; jrec < nrecs_; ++jrec) {
      for (std::size_t jvar = 0; jvar < vars_.size(); ++jvar) {
        const std::size_t jrecvar = jrec * vars_.size() + jvar;
        if (obs_num_[jrecvar] != missing_int) {
          for (std::size_t jp = 0; jp < prednames_.size(); ++jp) {
            const std::size_t index = this->index(jrec, jvar, jp);

            // Reset variances for bias predictor coeff. based on current data count
            if (obs_num_[jrecvar] <= minimal_required_obs_number_) {
              variances_[index] = smallest_variance_;
            }

            // Reset preconditioner L = \mathrm{A}^{-1}
            if (obs_num_[jrecvar] > 0)
              preconditioner_[index] = 1.0 / (1.0 / variances_[index] + ht_rinv_h_[index]);

            // Reset analysis variances
            if (obs_num_[jrecvar] > minimal_required_obs_number_) {
              if (ht_rinv_h_[index] > 0.0) {
                analysis_variances_[index] = 1.0 / (1.0 / variances_[index] + ht_rinv_h_[index]);
              } else {
                analysis_variances_[index] = largest_analysis_variance_;
              }
            }
          }
        } else {
          // Set variances, preconditioner and analysis variances to missing
          // for variables opted out of bias correction
          for (std::size_t jp = 0; jp < prednames_.size(); ++jp) {
            const std::size_t index = this->index(jrec, jvar, jp);
            variances_[index] = missing;
            preconditioner_[index] = missing;
            analysis_variances_[index] = missing;
          }
        }
      }
    }
  }
  oops::Log::trace() << "ObsBiasCovariance::linearize is done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasCovariance::multiply(const ObsBiasIncrement & dx1,
                                 ObsBiasIncrement & dx2) const {
  oops::Log::trace() << "ObsBiasCovariance::multiply starts" << std::endl;

  dx2.data().array() = dx1.data().array() * variances_.array();

  oops::Log::trace() << "ObsBiasCovariance::multiply is done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasCovariance::inverseMultiply(const ObsBiasIncrement & dx1,
                                        ObsBiasIncrement & dx2) const {
  oops::Log::trace() << "ObsBiasCovariance::inverseMultiply starts" << std::endl;

  dx2.data().array() = dx1.data().array() / variances_.array();

  oops::Log::trace() << "ObsBiasCovariance::inverseMultiply is done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsBiasCovariance::randomize(ObsBiasIncrement & dx) const {
  oops::Log::trace() << "ObsBiasCovariance::randomize starts" << std::endl;
  if (dx) {
    const double missing = util::missingValue<double>();
    static util::NormalDistribution<double> dist(variances_.size());
    for (std::size_t jj = 0; jj < variances_.size(); ++jj) {
      if (variances_[jj] != missing) {
        dx.data()[jj] = dist[jj] * std::sqrt(variances_[jj]);
      } else {
        dx.data()[jj] = missing;
      }
    }
  }
  oops::Log::trace() << "ObsBiasCovariance::randomize is done" << std::endl;
}
// -----------------------------------------------------------------------------

std::unique_ptr<ObsBiasPreconditioner> ObsBiasCovariance::preconditioner() const {
    return std::make_unique<ObsBiasPreconditioner> (preconditioner_);
}


// -----------------------------------------------------------------------------

}  // namespace ufo
