/*
 * (C) Copyright 2017-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/operators/crtm/ObsRadianceCRTM.h"

#include <algorithm>
#include <memory>
#include <utility>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Locations.h"
#include "oops/base/SamplingMethodSelector.h"
#include "oops/base/Variables.h"
#include "oops/interface/SampledLocations.h"
#include "oops/util/missingValues.h"

#include "ufo/fov/ReduceOverFieldOfView.h"
#include "ufo/fov/SampleFieldOfView.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/ObsTraits.h"
#include "ufo/operators/crtm/ObsRadianceCRTM.interface.h"
#include "ufo/SampledLocations.h"
#include "ufo/ScopedDefaultGeoVaLFormatChange.h"

namespace detail {

// If FOV-averaging is turned on, vars from this list will be averaged over the FOV ellipse
static const oops::Variables vars_to_fov_average(std::vector<std::string>{{
  "water_area_fraction",
  "land_area_fraction",
  "ice_area_fraction",
  "surface_snow_area_fraction",
  "surface_temperature_where_sea",
  "surface_temperature_where_land",
  "surface_temperature_where_ice",
  "surface_temperature_where_snow",
  "vegetation_area_fraction",
  "leaf_area_index",
  "volume_fraction_of_condensed_water_in_soil",
  "soil_temperature",
  "land_type_index_IGBP",
  "land_type_index_NPOESS",
  "land_type_index_USGS",
  "vegetation_type_index",
  "soil_type",
  "surface_snow_thickness",
  }});

// Selects method #0 for upper-air and non-FOV surface variables; #1 for surface FOV variables
class FovSelector : public oops::SamplingMethodSelector {
 public:
  FovSelector(const bool do_fov) : do_fov_average_(do_fov) {}

  size_t methodIndex(const std::string & varName) const override {
    if (do_fov_average_ && vars_to_fov_average.has(varName)) {
      // surface variables to average over FOV
      return 1;
    } else {
      // upper-air variables + surface variables with no averaging
      return 0;
    }
  }
 private:
  const bool do_fov_average_;
};

// This function gives the fallback values used to initialize CRTM data outside of JEDI masks, which
// could be either interpolation masks in GetValues or "area fraction" masks in FOV averaging.
//
// Empirically, the values used outside these masks don't appear to have any affect the CRTM
// radiance calculation, AS LONG AS they are reasonably physical values that don't cause the CRTM to
// crash. (So the simple alternative of keeping missingValues outside the GetValues masks doesn't
// work, because it leads to crashes inside CRTM.)
//
// So the actual values returned by this function, as long as reasonably physical, shouldn't matter.
// For now, we return the CRTM's default values for surface fields. The different values are
// collected into this single fallback function to help code readability.
double valueOutsideMask(const std::string & varname) {
  if (varname == "vegetation_type_index"
      || varname == "soil_type"
      || varname == "land_type_index_NPOESS"
      || varname == "land_type_index_IGBP") {
    return 1.0;
  } else if (varname == "leaf_area_index") {
    return 3.5;
  } else if (varname == "soil_temperature"
             || varname == "surface_temperature_where_sea"
             || varname == "surface_temperature_where_land"
             || varname == "surface_temperature_where_ice"
             || varname == "surface_temperature_where_snow") {
    return 283.0;  // K
  } else if (varname == "vegetation_area_fraction") {
    return 0.3;
  } else if (varname == "surface_snow_thickness") {
    return 50.0;  // mm
  } else if (varname == "volume_fraction_of_condensed_water_in_soil") {
    return 0.05;  // g/cm3
  } else {
    ABORT("Function valueOutsideMask doesn't yet provide a fallback for variable: " + varname);
    return util::missingValue<double>();
  }
}

}  // namespace detail

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsRadianceCRTM> makerCRTM_("CRTM");

// -----------------------------------------------------------------------------

ObsRadianceCRTM::ObsRadianceCRTM(const ioda::ObsSpace & odb,
                                 const Parameters_ & parameters)
  : ObsOperatorBase(odb), keyOperRadianceCRTM_(0),
    odb_(odb), varin_(), parameters_(parameters),
    do_fov_average_(parameters.DoFovAverage),
    fov_sample_resol_(parameters.FovSampleResol),
    sample_weights_()
{
  // parse channels from the config and create variable names
  const oops::Variables & observed = odb.assimvariables();
  std::vector<int> channels_list = observed.channels();

  // call Fortran setup routine
  ufo_radiancecrtm_setup_f90(keyOperRadianceCRTM_, parameters_.toConfiguration(),
                             channels_list.size(), channels_list[0], varin_, odb.comm());

  oops::Log::info() << "ObsRadianceCRTM channels: " << channels_list << std::endl;
  oops::Log::trace() << "ObsRadianceCRTM created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsRadianceCRTM::~ObsRadianceCRTM() {
  ufo_radiancecrtm_delete_f90(keyOperRadianceCRTM_);
  oops::Log::trace() << "ObsRadianceCRTM destructed" << std::endl;
}

// -----------------------------------------------------------------------------

// If FOV averaging is NOT requested, returns a single-element vector of obs locations.
// If FOV averaging is requested, returns a two-element vector:
// - the first element is the obs locations; this is intended to be paired with upper-air fields
//   (and any surface fields exempt from FOV averaging)
// - the second element is an expanded set of locations that sample the field of view ellipse
//   around each observation; this is intended to be paired with surface fields that should
//   receive FOV averaging. The resolution of the ellipse sampling is controlled via the yaml
//   option "fov sample points per semimajor axis": this generates a latlon grid around the obs,
//   from which points falling within the ellipse are kept.
ObsRadianceCRTM::Locations_ ObsRadianceCRTM::locations() const {
  oops::Log::trace() << "ObsRadianceCRTM locations started" << std::endl;

  typedef oops::SampledLocations<ObsTraits> SampledLocations_;

  std::vector<SampledLocations_> sampledLocations{};
  std::unique_ptr<oops::SamplingMethodSelector> selector{};

  // FOV center coords
  const size_t nlocs = odb_.nlocs();
  // floats are needed for the Locations constructor :(
  std::vector<float> lons(nlocs);
  std::vector<float> lats(nlocs);
  std::vector<util::DateTime> times(nlocs);
  odb_.get_db("MetaData", "longitude", lons);
  odb_.get_db("MetaData", "latitude", lats);
  odb_.get_db("MetaData", "dateTime", times);

  sampledLocations.push_back(SampledLocations_(
        std::make_unique<SampledLocations>(
          lons, lats, times, odb_.distribution())));

  if (do_fov_average_) {
    // Get attributes needed for GSI calls
    const ioda::ObsGroup og = odb_.getObsGroup();
    ASSERT(og.atts.exists("sensor"));
    ASSERT(og.atts.exists("platform"));
    const std::string & sensor = og.atts.read<std::string>("sensor");
    const std::string & platform = og.atts.read<std::string>("platform");

    // Get observation data that will defined the field-of-view
    std::vector<int> scan_positions(nlocs);
    std::vector<float> azimuth_angles(nlocs);
    odb_.get_db("MetaData", "sensorScanPosition", scan_positions);
    odb_.get_db("MetaData", "sensorAzimuthAngle", azimuth_angles);

    // Sample the FOV ellipse
    std::vector<float> sample_lons;
    std::vector<float> sample_lats;
    std::vector<util::DateTime> sample_times;
    std::vector<util::Range<size_t>> sample_ranges;  // range [begin,end) of samples for each obs
    fov::getSampleLocationsAndWeights(
        sample_lons, sample_lats, sample_times, sample_ranges, sample_weights_,
        sensor, platform, fov_sample_resol_,
        lons, lats, times, scan_positions, azimuth_angles);

    sampledLocations.push_back(SampledLocations_(
          std::make_unique<SampledLocations>(
            sample_lons, sample_lats, sample_times,
            odb_.distribution(), sample_ranges)));
  }

  oops::Log::trace() << "ObsRadianceCRTM locations done" << std::endl;
  return Locations_(std::move(sampledLocations),
                    std::make_unique<detail::FovSelector>(do_fov_average_));
}

// -----------------------------------------------------------------------------

void ObsRadianceCRTM::computeReducedVars(const oops::Variables & vars, GeoVaLs & geovals) const {
  oops::Log::trace() << "ObsRadianceCRTM computeReducedVars started" << std::endl;

  // Sanity check: GeoVaLs do not yet contain reduced vars
  // (OR: reduced vars match sampled vars; this is not production behavior, but needed for tests
  //      because geovals files are read in as sampled AND reduced variables...)
  ASSERT(geovals.getReducedVars().size() == 0 || geovals.getReducedVars() == geovals.getVars());

  // Sanity check:
  // - all vars to reduce are present in geovals
  // - any vars to reduce that were not requested by CRTM are not FOV vars (simplifies logic later)
  const oops::Variables & gvars = geovals.getVars();
  for (const auto & v : vars.variables()) {
    ASSERT(gvars.has(v));
  }
  oops::Variables vars_not_for_crtm = vars;
  vars_not_for_crtm -= varin_;
  for (const auto & v : vars_not_for_crtm.variables()) {
    ASSERT(!detail::vars_to_fov_average.has(v));
  }

  // Allocate reduced vars in geovals
  // (UNLESS: reduced vars match sampled vars, as can occur in tests; then do nothing)
  std::vector<size_t> sizes(gvars.size());
  for (size_t i = 0; i < gvars.size(); ++i) {
    sizes[i] = geovals.nlevs(gvars[i], GeoVaLFormat::SAMPLED);
  }
  if (geovals.getReducedVars().size() == 0) {
    geovals.addReducedVars(gvars, sizes);
  }

  // Sanity check: non-FOV vars are aliased
  for (size_t i = 0; i < gvars.size(); ++i) {
    if (!detail::vars_to_fov_average.has(gvars[i])) {  // non-FOV var
      if (!geovals.areReducedAndSampledFormatsAliased(gvars[i])) {
        ABORT("Internal logic error: expected aliased reduced variable " + gvars[i]);
      }
    }
  }

  // Fill the reduced vars...
  if (do_fov_average_) {
    // ... by averaging each obs over its FOV samples
    fillReducedVarsByMaskedAveraging(geovals);
  } else {
    // ... by looking for missingValues in the geovals and masking these out with a fallback value.
    // This branch is to handle interpolation masks in GetValues, which can lead to missingValues
    // for GeoVaL variables falling outside of their interpolation mask. This branch is doing work
    // that's not strictly needed when no masks are enabled, but the work should not be too costly.
    fillReducedVarsByMaskedCopy(geovals);
  }

  oops::Log::trace() << "ObsRadianceCRTM computeReducedVars done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceCRTM::fillReducedVarsByMaskedAveraging(GeoVaLs & geovals) const {
  oops::Log::trace() << "ObsRadianceCRTM fillReducedVarsByMaskedAveraging started" << std::endl;

  oops::Variables vars_to_reduce = varin_;
  vars_to_reduce.intersection(detail::vars_to_fov_average);

  // Sample ranges to use when reducing each geoval
  std::vector<util::Range<size_t>> sample_ranges;
  geovals.getProfileIndicesGroupedByLocation(vars_to_reduce[0], sample_ranges);
  const size_t nlocs = sample_ranges.size();
  // The statement below assumes the samples run over obs locs in increasing order
  const size_t nsamples = sample_ranges.back().end;

  const size_t expected_field_counter = vars_to_reduce.size();
  size_t field_counter = 0;

  // Helper to compute area fractions (averaged per obs) and masks (per sample) associated with
  // a particular surface type.
  const auto & maskHelper = [this, &nsamples, &nlocs, &sample_ranges,
                             &geovals, &field_counter](
      const std::string & varname,
      std::vector<double> & sample_values,
      std::vector<double> & average) -> void {
    ASSERT(geovals.nprofiles(varname, GeoVaLFormat::SAMPLED) == nsamples);
    ASSERT(sample_values.size() == nsamples);
    ASSERT(average.size() == nlocs);
    geovals.getAtLevel(sample_values, varname, 0, GeoVaLFormat::SAMPLED);
    fov::average(average, sample_ranges, sample_values, sample_weights_);
    geovals.putAtLevel(average, varname, 0, GeoVaLFormat::REDUCED);
    field_counter++;
  };

  std::vector<double> water_mask(nsamples);
  std::vector<double> water_avg(nlocs);
  maskHelper("water_area_fraction", water_mask, water_avg);

  std::vector<double> land_mask(nsamples);
  std::vector<double> land_avg(nlocs);
  maskHelper("land_area_fraction", land_mask, land_avg);

  std::vector<double> ice_mask(nsamples);
  std::vector<double> ice_avg(nlocs);
  maskHelper("ice_area_fraction", ice_mask, ice_avg);

  std::vector<double> snow_mask(nsamples);
  std::vector<double> snow_avg(nlocs);
  maskHelper("surface_snow_area_fraction", snow_mask, snow_avg);

  // Sanity check the area fractions add to 1 for each obs
  for (size_t i = 0; i < nlocs; ++i) {
    const double check = water_avg[i] + land_avg[i] + ice_avg[i] + snow_avg[i];
    if (std::fabs(check - 1.0) > 1e-5) {  // arbitrary threshold
      ABORT("Invalid surface area fractions when computing FOV averages for CRTM radiance");
    }
  }

  // Helper to average over FOV using a mask
  const auto & averageHelper = [this, &nsamples, &nlocs, &sample_ranges,
                                &geovals, &field_counter](
      const std::string & varname,
      const std::vector<double> mask,
      const std::vector<double> sample_mask) -> void {
    std::vector<double> sample_values(nsamples);
    std::vector<double> average(nlocs);
    geovals.getAtLevel(sample_values, varname, 0, GeoVaLFormat::SAMPLED);
    const double valueOutsideMask = detail::valueOutsideMask(varname);
    fov::averageWithMask(average, sample_ranges, sample_values, sample_weights_,
                         mask, sample_mask, valueOutsideMask);
    geovals.putAtLevel(average, varname, 0, GeoVaLFormat::REDUCED);
    field_counter++;
  };

  averageHelper("surface_temperature_where_sea", water_avg, water_mask);
  averageHelper("surface_temperature_where_land", land_avg, land_mask);
  averageHelper("surface_temperature_where_ice", ice_avg, ice_mask);
  averageHelper("surface_temperature_where_snow", snow_avg, snow_mask);

  averageHelper("leaf_area_index", land_avg, land_mask);
  averageHelper("soil_temperature", land_avg, land_mask);
  averageHelper("vegetation_area_fraction", land_avg, land_mask);
  averageHelper("volume_fraction_of_condensed_water_in_soil", land_avg, land_mask);
  averageHelper("surface_snow_thickness", snow_avg, snow_mask);

  // Helper to identify the dominant surface classification type using a mask
  const auto & surfaceTypeHelper = [this, &nsamples, &nlocs, &sample_ranges,
                                    &geovals, &field_counter](
      const std::string & varname,
      const std::vector<double> mask,
      const std::vector<double> sample_mask) -> void {
    std::vector<double> sample_int_values(nsamples);
    std::vector<double> dominant_int(nlocs);
    geovals.getAtLevel(sample_int_values, varname, 0, GeoVaLFormat::SAMPLED);
    const int valueOutsideMask = static_cast<int>(detail::valueOutsideMask(varname));
    fov::dominantIntegerValue(dominant_int, sample_ranges, sample_int_values, sample_weights_,
                              mask, sample_mask, valueOutsideMask);
    geovals.putAtLevel(dominant_int, varname, 0, GeoVaLFormat::REDUCED);
    field_counter++;
  };

  if (varin_.has("vegetation_type_index") && varin_.has("soil_type")) {
    surfaceTypeHelper("vegetation_type_index", land_avg, land_mask);
    surfaceTypeHelper("soil_type", land_avg, land_mask);
  } else if (varin_.has("land_type_index_NPOESS")) {
    surfaceTypeHelper("land_type_index_NPOESS", land_avg, land_mask);
  } else if (varin_.has("land_type_index_IGBP")) {
    surfaceTypeHelper("land_type_index_IGBP", land_avg, land_mask);
  } else {
    ABORT("Inconsistent or unsupported surface types");
  }

  // Check every FOV variable was handled
  ASSERT(field_counter == expected_field_counter);
  oops::Log::trace() << "ObsRadianceCRTM fillReducedVarsByMaskedAveraging done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceCRTM::fillReducedVarsByMaskedCopy(GeoVaLs & geovals) const {
  oops::Log::trace() << "ObsRadianceCRTM fillReducedVarsByMaskedCopy started" << std::endl;

  // Helper to replace missing values with a fallback value
  // (This could be made more robust by checking missing values only occur when the "correct" area
  //  fraction is zero: e.g., soil_type should only be missing if land_area_fraction is zero.
  //  Ideally this would be done in debug builds only, so requires some prior infrastructure.)
  const auto & maskHelper = [&geovals](const std::string & name) -> void {
    const size_t nlocs = geovals.nlocs();
    std::vector<double> vals(nlocs);
    geovals.getAtLevel(vals, name, 0, GeoVaLFormat::SAMPLED);
    bool fixed = false;
    for (size_t i = 0; i < nlocs; ++i) {
      if (vals[i] == util::missingValue<double>()) {
        vals[i] = detail::valueOutsideMask(name);
        fixed = true;
      }
    }
    if (fixed) geovals.putAtLevel(vals, name, 0, GeoVaLFormat::REDUCED);
  };

  maskHelper("surface_temperature_where_sea");
  maskHelper("surface_temperature_where_land");
  maskHelper("surface_temperature_where_ice");
  maskHelper("surface_temperature_where_snow");

  maskHelper("leaf_area_index");
  maskHelper("soil_temperature");
  maskHelper("vegetation_area_fraction");
  maskHelper("volume_fraction_of_condensed_water_in_soil");
  maskHelper("surface_snow_thickness");

  if (varin_.has("vegetation_type_index") && varin_.has("soil_type")) {
    maskHelper("vegetation_type_index");
    maskHelper("soil_type");
  } else if (varin_.has("land_type_index_NPOESS")) {
    maskHelper("land_type_index_NPOESS");
  } else if (varin_.has("land_type_index_IGBP")) {
    maskHelper("land_type_index_IGBP");
  } else {
    ABORT("Inconsistent or unsupported surface types");
  }

  oops::Log::trace() << "ObsRadianceCRTM fillReducedVarsByMaskedCopy done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceCRTM::simulateObs(const GeoVaLs & gom, ioda::ObsVector & ovec,
                                  ObsDiagnostics & dvec, const QCFlags_t & qc_flags) const {
  oops::Log::trace() << "ObsRadianceCRTM simulateObs started" << std::endl;

  // Write to geovals file after calling computeReducedVars; this can be useful to check
  // the outputs of the reduction before it passes through the complicated CRTM.
  if (parameters_.gvOut.value() != boost::none) {
    const eckit::LocalConfiguration & config = parameters_.gvOut.value().value();
    gom.write(config);
  }

  // Always pass the reduced format to the CRTM; even when not doing FOV averaging, this allows
  // for fixing of the GeoVaLs to accomodate interpolation masks in GetValues.
  ScopedDefaultGeoVaLFormatChange changeFormat(gom, GeoVaLFormat::REDUCED);

  ufo_radiancecrtm_simobs_f90(keyOperRadianceCRTM_, gom.toFortran(), odb_,
                              ovec.nvars(), ovec.nlocs(), ovec.toFortran(),
                              dvec.toFortran(),
                              reinterpret_cast<const void*>(&qc_flags));

  oops::Log::trace() << "ObsRadianceCRTM simulateObs done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsRadianceCRTM::print(std::ostream & os) const {
  os << "ObsRadianceCRTM::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
