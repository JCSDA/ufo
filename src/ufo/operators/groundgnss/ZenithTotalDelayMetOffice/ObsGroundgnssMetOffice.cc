/*
 * (C) Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/groundgnss/ZenithTotalDelayMetOffice/ObsGroundgnssMetOffice.h"

#include <memory>
#include <numeric>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Locations.h"
#include "oops/base/SamplingMethodSelector.h"
#include "oops/base/Variables.h"
#include "oops/interface/SampledLocations.h"
#include "oops/util/Logger.h"
#include "oops/util/Range.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/ObsTraits.h"
#include "ufo/SampledLocations.h"

namespace ufo {

namespace {

// Selects the FOV sampling method (#1) for the variables this operator reads and the single-path
// method (#0) for everything else. ObsDiagnostics (e.g. bias-correction terms) and filter or
// bias-predictor GeoVaLs are then sized by nlocs, as they expect, rather than by nlocs * num_fov.
class FovSelector : public oops::SamplingMethodSelector {
 public:
  explicit FovSelector(const oops::Variables & fovVars) : fovVars_(fovVars) {}

  size_t methodIndex(const oops::Variable & var) const override {
    return fovVars_.has(var) ? 1 : 0;
  }

 private:
  const oops::Variables fovVars_;
};

}  // namespace

// ----------------------------------------------------------------------------
static ObsOperatorMaker<ObsGroundgnssMetOffice> makerGroundgnssMetOffice_("GroundgnssMetOffice");
// -----------------------------------------------------------------------------

ObsGroundgnssMetOffice::ObsGroundgnssMetOffice(const ioda::ObsSpace & odb,
                                       const Parameters_ & parameters)
  : ObsOperatorBase(odb), odb_(odb), numfov_(parameters.numFOV.value())
{
  const std::vector<std::string> vv{"air_pressure_levels",
                                    "water_vapor_mixing_ratio_wrt_moist_air",
                                    "geopotential_height", "geopotential_height_levels"};
  varin_.reset(new oops::Variables(vv));

  ufo_groundgnss_metoffice_setup_f90(keyOperGroundgnssMetOffice_,
                                     parameters.vertInterpOPS.value(),
                                     parameters.pseudoLevels.value(),
                                     parameters.minTempGrad.value(),
                                     parameters.dryRefractivityConstant.value(),
                                     parameters.wetRefractivityConstant.value(),
                                     parameters.numFOV.value(),
                                     parameters.fovRadius.value());

  oops::Log::trace() << "ObsGroundgnssMetOffice constructor done." << std::endl;
}

// -----------------------------------------------------------------------------

ObsGroundgnssMetOffice::~ObsGroundgnssMetOffice() {
  ufo_groundgnss_metoffice_delete_f90(keyOperGroundgnssMetOffice_);
  oops::Log::trace() << "ObsGroundgnssMetOffice destructor done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGroundgnssMetOffice::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                                     ObsDiagnostics & d, const QCFlags_t & qc_flags) const {
  ufo_groundgnss_metoffice_simobs_f90(keyOperGroundgnssMetOffice_, gom.toFortran(), odb_,
                                  ovec.size(), ovec.toFortran());
  oops::Log::trace() << "ObsGroundgnssMetOffice::simulateObs done." << std::endl;
}

// -----------------------------------------------------------------------------

ObsGroundgnssMetOffice::Locations_
ObsGroundgnssMetOffice::locations() const {
  using SampledLocations_ = oops::SampledLocations<ObsTraits>;

  // Without FOV averaging each location is sampled along a single path, which is exactly what
  // the base class already does.
  if (numfov_ <= 1) {
    return ObsOperatorBase::locations();
  }

  const size_t nlocs = odb_.nlocs();
  const size_t npaths = nlocs * numfov_;

  // Method #0: one path per location, at the observation itself, as in the base class.
  std::vector<float> obsLons(nlocs);
  std::vector<float> obsLats(nlocs);
  std::vector<util::DateTime> obsTimes(nlocs);
  odb_.get_db("MetaData", "longitude", obsLons);
  odb_.get_db("MetaData", "latitude", obsLats);
  odb_.get_db("MetaData", "dateTime", obsTimes);

  std::vector<SampledLocations_> sampledLocations;
  sampledLocations.push_back(SampledLocations_(
    std::make_unique<SampledLocations>(obsLons, obsLats, obsTimes, odb_.distribution())));

  // Method #1: num_fov paths per location spread over the field of view.
  std::vector<float> lons(npaths);
  std::vector<float> lats(npaths);
  std::vector<util::DateTime> times(npaths);

  // Every sampling point of an observation inherits that observation's time.
  for (size_t jloc = 0; jloc < nlocs; ++jloc) {
    for (size_t jfov = 0; jfov < numfov_; ++jfov) {
      times[jloc*numfov_ + jfov] = obsTimes[jloc];
    }
  }

  // The sampling pattern itself is laid out in Fortran, next to the code that consumes it.
  if (npaths > 0) {
     ufo_groundgnss_metoffice_locs_init_f90(keyOperGroundgnssMetOffice_, odb_, npaths,
                                         lons[0], lats[0]);
  }

  std::vector<util::Range<size_t>> pathsGroupedByLocation(nlocs);
  for (size_t jloc = 0; jloc < nlocs; ++jloc) {
    pathsGroupedByLocation[jloc].begin = jloc * numfov_;
    pathsGroupedByLocation[jloc].end = (jloc + 1) * numfov_;
  }

  sampledLocations.push_back(SampledLocations_(
    std::make_unique<SampledLocations>(lons, lats, times, odb_.distribution(),
                                       std::move(pathsGroupedByLocation))));

  return Locations_(std::move(sampledLocations), std::make_unique<FovSelector>(*varin_));
}

// -----------------------------------------------------------------------------

void ObsGroundgnssMetOffice::computeReducedVars(const oops::Variables & vars,
                                                GeoVaLs & geovals) const {
  // Without FOV averaging the sampled and reduced formats are aliased, so the base class
  // handles this correctly.
  if (numfov_ <= 1) {
    ObsOperatorBase::computeReducedVars(vars, geovals);
    return;
  }

  if (odb_.nlocs() == 0) {
    return;
  }

  // Add a reduced format for every sampled variable, not just those in `vars`: filters and bias
  // predictors do not always request explicitly all the variables they go on to read.
  const oops::Variables & gvars = geovals.getVars();
  oops::Variables missingVars = gvars;
  missingVars -= geovals.getReducedVars();
  if (missingVars.size() != 0) {
    std::vector<size_t> missingSizes(missingVars.size());
    for (size_t ivar = 0; ivar < missingVars.size(); ++ivar) {
      missingSizes[ivar] = geovals.nlevs(missingVars[ivar], GeoVaLFormat::SAMPLED);
    }
    geovals.addReducedVars(missingVars, missingSizes);
  }

  // Average every variable sampled along several paths per location. Variables on the single-path
  // method share storage between the two formats, so they need nothing further.
  for (size_t ivar = 0; ivar < gvars.size(); ++ivar) {
    const oops::Variable & var = gvars[ivar];
    if (geovals.areReducedAndSampledFormatsAliased(var)) continue;

    std::vector<util::Range<size_t>> sample_ranges;
    geovals.getProfileIndicesGroupedByLocation(var, sample_ranges, GeoVaLFormat::SAMPLED);
    const size_t nlocs = sample_ranges.size();
    std::vector<double> sample_values(geovals.nprofiles(var, GeoVaLFormat::SAMPLED));
    std::vector<double> average(nlocs);

    for (size_t jlev = 0; jlev < geovals.nlevs(var, GeoVaLFormat::SAMPLED); ++jlev) {
      geovals.getAtLevel(sample_values, var, jlev, GeoVaLFormat::SAMPLED);
      for (size_t jloc = 0; jloc < nlocs; ++jloc) {
        // Equal weighting of the FOV sampling points, matching the Fortran forward operator.
        const double sum = std::accumulate(sample_values.begin() + sample_ranges[jloc].begin,
                                           sample_values.begin() + sample_ranges[jloc].end, 0.0);
        const size_t count = sample_ranges[jloc].end - sample_ranges[jloc].begin;
        average[jloc] = sum / static_cast<double>(count);
      }
      geovals.putAtLevel(average, var, jlev, GeoVaLFormat::REDUCED);
    }
  }

  oops::Log::trace() << "ObsGroundgnssMetOffice::computeReducedVars done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsGroundgnssMetOffice::print(std::ostream & os) const {
  os << "ObsGroundgnssMetOffice::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
