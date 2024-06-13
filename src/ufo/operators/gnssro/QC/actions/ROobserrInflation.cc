/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/gnssro/QC/actions/ROobserrInflation.h"
#include <algorithm>
#include <set>

#include "ioda/ObsDataVector.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/utils/StringUtils.h"

namespace ufo {

// -----------------------------------------------------------------------------

static FilterActionMaker<ROobserrInflation> makerInflateErr_("RONBAMErrInflate");

// -----------------------------------------------------------------------------

ROobserrInflation::ROobserrInflation(const Parameters_ &)
  : allvars_() {
}
// -----------------------------------------------------------------------------

void ROobserrInflation::apply(const Variables & vars,
                              const std::vector<std::vector<bool>> & flagged,
                              ObsFilterData & data,
                              int filterQCflag,
                              ioda::ObsDataVector<int> & flags,
                              ioda::ObsDataVector<float> & obserr) const {
  const float missing = util::missingValue<float>();
  size_t nlocs = data.nlocs();
  ioda::ObsDataVector<int> layeridx(data.obsspace(), "modelLayerIndex", "ObsDiag");
  ioda::ObsDataVector<int> recordidx(data.obsspace(), "RecordNumberIndex", "ObsDiag");
  const std::vector<size_t> & recordNumbers = data.obsspace().recnum();

  std::vector<int> rec_idx(nlocs);
  std::vector<int> layer_idx(nlocs);
  std::vector<float> factor(nlocs);
  const int maxlev = 500;
  std::vector<std::vector<int> > super_obs_inlayer(maxlev, std::vector<int>(nlocs));

  for (size_t i = 0; i < maxlev; ++i) {
      for (size_t j = 0; j < nlocs; ++j) {
          super_obs_inlayer[i][j] = 0;
      }
  }
  int irec = 1;
  for (size_t jobs = 0; jobs < nlocs; ++jobs) {
      if (jobs > 0 && recordNumbers[jobs] != recordNumbers[jobs-1]) irec = irec +1;
      rec_idx[jobs] = recordidx[0][jobs] - 1;
      layer_idx[jobs] = layeridx[0][jobs];

      if (flags[0][jobs] == 0)  super_obs_inlayer[layer_idx[jobs]][rec_idx[jobs]]++;
  }
  for (size_t jobs = 0; jobs < nlocs; ++jobs) {
      factor[jobs] = 1.0;
      if (super_obs_inlayer[layer_idx[jobs]][rec_idx[jobs]] > 0) {
         factor[jobs] = super_obs_inlayer[layer_idx[jobs]][rec_idx[jobs]];
         factor[jobs] = sqrt(factor[jobs]);
      }
      if (obserr[0][jobs] != missing) obserr[0][jobs] *= factor[jobs];
  }
}
// -----------------------------------------------------------------------------

}  // namespace ufo
