/*
 * (C) Copyright 2025 Space Sciences and Engineering, LLC (dba PlanetiQ).
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * author Steve Marshall (smarshall@planetiq.com)
 */
#include "oops/util/Logger.h"
#include "ufo/operators/gnssro/utils/GnssroProfileExtractor.h"

namespace ufo {

// -----------------------------------------------------------------------------

GnssroProfileExtractor::GnssroProfileExtractor(const ioda::ObsSpace & odb)
  : slices_()
{
  oops::Log::trace() << "GnssroProfileExtractor() created" << std::endl;
  std::size_t sumSliceSizes = 0;
  if (odb.nlocs() > 0)
  {
    // Get sequenceNumber from the observations.
    std::vector<int> seqNums(odb.nlocs());
    std::vector<util::DateTime> times(odb.nlocs());
    odb.get_db("MetaData", "sequenceNumber", seqNums);
    odb.get_db("MetaData", "dateTime", times);

    //  If the data is empty, there is nothing to do.
    if (seqNums.empty()) {
      oops::Log::trace() << "GnssroProfileExtractor(): constructor complete: "
                         << "found no sequenceNumber values" << std::endl;
      return;
    }
    oops::Log::debug() << "GnssroProfileExtractor(): found " << seqNums.size()
                       << " sequenceNumber values" << std::endl;

    // Find the range associated with each unique sequence number.
    // NOTE: due to filtering and distribution of work between ranks,
    // there may be gaps in the integer values and the sequence numbers
    // will not be in order.
    size_t start = 0;
    int seqNum = seqNums[0];
    util::DateTime profileTime = times[0];
    for (std::size_t iloc = 1; iloc < seqNums.size(); ++iloc) {
      if (seqNums[iloc] != seqNum) {
        slices_.push_back(GnssroProfileSlice(seqNum, profileTime, start, iloc));
        start = iloc;
        seqNum = seqNums[iloc];
        profileTime = times[iloc];
      }
    }
    slices_.push_back(GnssroProfileSlice(seqNum, profileTime, start, seqNums.size()));
    sumSliceSizes += seqNums.size();
  }

  oops::Log::trace() << "GnssroProfileExtractor() constructor complete: " << slices_.size()
                     << " profiles found with total of " << sumSliceSizes << " obs; nlocs="
                     << odb.nlocs() << " expected" << std::endl;
}
// -----------------------------------------------------------------------------

GnssroProfileExtractor::~GnssroProfileExtractor()
{
}
// -----------------------------------------------------------------------------

}  // namespace ufo

