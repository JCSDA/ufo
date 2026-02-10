/*
 * (C) Copyright 2023 United States Government as represented by the Administrator of the
 * National Aeronautics and Space Administration. All Rights Reserved.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/ObsErrorFactorDuplicateCheck.h"

#include <float.h>
#include <math.h>

#include <algorithm>
#include <cmath>

#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"
#include "ioda/ObsDataVector.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsAccessor.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/GeoVaLs.h"
#include "ufo/utils/Constants.h"
#include "ufo/utils/PiecewiseLinearInterpolation.h"
#include "ufo/variabletransforms/Formulas.h"

namespace ufo {

static ObsFunctionMaker<ObsErrorFactorDuplicateCheck> makerSteps_("ObsErrorFactorDuplicateCheck");

// -----------------------------------------------------------------------------

ObsErrorFactorDuplicateCheck::ObsErrorFactorDuplicateCheck(const eckit::Configuration &config)
  : invars_() {
  oops::Log::trace() << "ObsErrorFactorDuplicateCheck constructor" << std::endl;
  // Initialize options
  options_.reset(new ObsErrorFactorDuplicateCheckParameters());
  options_->deserialize(config);
  const std::string errgrp = options_->original_obserr.value();
  const std::string flaggrp = options_->testQCflag.value();
  const bool use_air_pres = options_->use_air_pressure.value();
  const std::string var = options_->variable.value();

  // Include list of required data from MetaData
  invars_ += Variable(errgrp+"/"+var);
  invars_ += Variable(flaggrp+"/"+var);
  invars_ += Variable("MetaData/latitude");
  invars_ += Variable("MetaData/longitude");
  invars_ += Variable("MetaData/dateTime");
  invars_ += Variable("MetaData/pressure");
}

// -----------------------------------------------------------------------------

ObsErrorFactorDuplicateCheck::~ObsErrorFactorDuplicateCheck() {
  oops::Log::trace() << "ObsErrorFactorDuplicateCheck destructor" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename T>
std::vector<T> ObsErrorFactorDuplicateCheck::getGlobalVariable(const ObsFilterData & data,
                                                      const std::string &group,
                                                      const std::string &variable) const {
  auto & obsdb = data.obsspace();
  std::vector<T> values(obsdb.nlocs());
  data.get(Variable(group+"/"+variable), values);
  obsdb.distribution()->allGatherv(values);
  return values;
}

int ObsErrorFactorDuplicateCheck::getDispl(const ObsFilterData & data,
                                           const std::vector<int> & rank_nlocs) const {
  // Retrieves displacement for a given rank
  auto & obsdb = data.obsspace();
  int rank = obsdb.comm().rank();
  if (rank == 0) return 0;
  return (accumulate(rank_nlocs.begin(), (rank_nlocs.begin()+(rank)), 0));
}

// -----------------------------------------------------------------------------
void ObsErrorFactorDuplicateCheck::compute(const ObsFilterData & data,
                                     ioda::ObsDataVector<float> & obserr) const {
  oops::Log::trace() << "ObsErrorFactorDuplicateCheck compute start" << std::endl;
  auto & obsdb = data.obsspace();
  int nlocs = data.nlocs();
  int commsize = obsdb.comm().size();
  std::vector<int> root_rank_nlocs(commsize);
  std::vector<int> rank_nlocs(commsize);

  const std::string var = options_->variable.value();
  const bool use_air_pres = options_->use_air_pressure.value();

  // Get displacements
  obsdb.comm().allGather(nlocs, rank_nlocs.begin(), rank_nlocs.end());
  int displ = getDispl(data, rank_nlocs);

  // Get local and global arrays for variables
  std::vector<float> lat_local(nlocs);
  data.get(Variable("MetaData/latitude"), lat_local);
  std::vector<float> lat_global = getGlobalVariable<float>(data, "MetaData", "latitude");

  std::vector<float> lon_local(nlocs);
  data.get(Variable("MetaData/longitude"), lon_local);
  std::vector<float> lon_global = getGlobalVariable<float>(data, "MetaData", "longitude");

  std::vector<util::DateTime> datetime_local(nlocs);
  data.get(Variable("MetaData/dateTime"), datetime_local);
  std::vector<util::DateTime> datetime_global =
                          getGlobalVariable<util::DateTime>(data, "MetaData", "dateTime");

  std::vector<float> pres_local(nlocs);
  data.get(Variable("MetaData/pressure"), pres_local);
  std::vector<float> pres_global = getGlobalVariable<float>(data, "MetaData", "pressure");

  std::vector<int> qcflag_local(nlocs);
  const std::string flaggrp = options_->testQCflag.value();
  data.get(Variable(flaggrp+"/"+var), qcflag_local);
  std::vector<int> qcflag_global = getGlobalVariable<int>(data, flaggrp, var);

  std::vector<float> currentObserr(nlocs);
  const std::string errgrp = options_->original_obserr.value();
  data.get(Variable(errgrp+"/"+var), currentObserr);

  // Duplicate set up
  std::vector<double> dup(qcflag_local.size(), 1);
  std::vector<int> inds;
  std::vector<int> inds_global;
  double tfact;
  double duplicate_retention_weight = options_->duplicate_retention_weight.value();
  double time_window_hours = options_->time_window_hours.value();
  double one = 1.0;
  double hours_to_seconds = 3600.0;

  // Find indicies that pass through qc for global and local arrs
  for (size_t i = 0; i < nlocs; ++i) {
    if (qcflag_local[i] == 0) inds.push_back(i);
    obserr[0][i] = 1.0;
  }
  for (size_t i = 0; i < qcflag_global.size(); ++i) {
    if (qcflag_global[i] == 0) inds_global.push_back(i);
  }

  int qc_len_local = inds.size();
  int qc_len_global = inds_global.size();
  int obserr_rank = obsdb.comm().rank();

  // Option to use pressure
  if (use_air_pres) {
    for (size_t iobs = 0; iobs < qc_len_local; ++iobs) {
      for (size_t jobs = 0; jobs < qc_len_global; ++jobs) {
        if ( (lat_local[inds[iobs]] == lat_global[inds_global[jobs]]) &&
           (lon_local[inds[iobs]] == lon_global[inds_global[jobs]]) &&
           (pres_local[inds[iobs]] == pres_global[inds_global[jobs]])) {
            // If not equal to itself
            if (inds_global[jobs] != (inds[iobs] + displ)) {
               tfact = std::min(one, (abs(((datetime_local[inds[iobs]] -
                                      datetime_global[inds_global[jobs]]).toSeconds()))/
                                      hours_to_seconds)/time_window_hours);
               dup[inds[iobs]] = dup[inds[iobs]]+1.0-tfact*tfact*(1.0-duplicate_retention_weight);
            }
        }
      }
      if (dup[inds[iobs]] > 1 ) obserr[0][inds[iobs]] = sqrt(dup[inds[iobs]]);
  }
  } else {
    for (size_t iobs = 0; iobs < qc_len_local; ++iobs) {
      for (size_t jobs = 0; jobs < qc_len_global; ++jobs) {
        if ( (lat_local[inds[iobs]] == lat_global[inds_global[jobs]]) &&
             (lon_local[inds[iobs]] == lon_global[inds_global[jobs]]) ) {
            if (inds_global[jobs] != (inds[iobs] + displ)) {
               tfact = std::min(one, (abs(((datetime_local[inds[iobs]] -
                                      datetime_global[inds_global[jobs]]).toSeconds()))/
                                      hours_to_seconds)/time_window_hours);
               dup[inds[iobs]] = dup[inds[iobs]]+1.0-tfact*tfact*(1.0-duplicate_retention_weight);
            }
        }
      }
      if ( dup[inds[iobs]] > 1 ) obserr[0][inds[iobs]] = sqrt(dup[inds[iobs]]);
    }
  }
  oops::Log::trace() << "ObsErrorFactorDuplicateCheck compute complete" << std::endl;
}
// -----------------------------------------------------------------------------

const ufo::Variables & ObsErrorFactorDuplicateCheck::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
