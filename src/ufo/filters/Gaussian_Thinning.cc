/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/Gaussian_Thinning.h"

#include <cmath>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "oops/util/Random.h"
#include "ufo/utils/Constants.h"

namespace ufo {

// -----------------------------------------------------------------------------

Gaussian_Thinning::Gaussian_Thinning(ioda::ObsSpace & obsdb, const eckit::Configuration & config,
                                     boost::shared_ptr<ioda::ObsDataVector<int> > flags,
                                     boost::shared_ptr<ioda::ObsDataVector<float> > obserr)
  : FilterBase(obsdb, config, flags, obserr)
{
  oops::Log::debug() << "Gaussian_Thinning: config = " << config_ << std::endl;
}

// -----------------------------------------------------------------------------

Gaussian_Thinning::~Gaussian_Thinning() {}

// -----------------------------------------------------------------------------

int Gaussian_Thinning::ll_to_idx(float in_lon, float in_lat, size_t n_lats,
                                 const std::vector<size_t> &n_lons) {
  // function to get i from lon/lat
  int ilat__ = 0;
  int ilon__ = 0;

  // ilat__ based on linear scaling from -90 to 90.  Equivalent to:
  // ilat__ = (in_lat - (-90.0) ) / (90.0 - (-90.0)) * n_lats;
  ilat__ = (in_lat - (-90.0) ) / (180.0) * n_lats;
  ilon__ = 0;
  ilon__ = std::accumulate(n_lons.begin(), n_lons.begin()+ilat__, 0);
  oops::Log::debug() << "Gaussian_Thinning: ilon__(acc): " << ilon__ << std::endl;

  // ilon__ based on linear scaling of lonmin=0, lonmax=360.  Equivalent to:
  // ilon__ = ilon__ + (in_lon - 0.0) / (360.0 - 0.0) * n_lons[ilat__];
  ilon__ += (in_lon) / (360.0) * n_lons[ilat__];
  return ilon__;
}

// -----------------------------------------------------------------------------

int Gaussian_Thinning::pres_to_idx(float in_pres, size_t n_vmesh,
                                   float vertical_mesh, float vertical_max) {
  // function to go from pres to k (vert. index) from pressure
  int idx__ = 0;
  if (n_vmesh > 1) {
     float new_vmin = vertical_max - n_vmesh * vertical_mesh;  // adjust vertical mesh for
                                                               //    unbalanced bin at top
     idx__ = (in_pres - vertical_max) / (new_vmin - vertical_max) * (n_vmesh);

     if (idx__ > n_vmesh - 1) { idx__ = n_vmesh - 1 ;}  // assign to top bin if above mesh
     if (idx__ < 0)           { idx__ = 0           ;}  // assign to bottom bin if below
  }
  return idx__;
}

// -----------------------------------------------------------------------------

int Gaussian_Thinning::dist_to_centroid(float ob_lon, float ob_lat, float c_lon, float c_lat) {
  // function to calculate distance from centroid
  const float deg2rad = Constants::deg2rad;
  const float re = Constants::mean_earth_rad;  // km

  float ob_lonr__ = ob_lon * deg2rad;
  float ob_latr__ = ob_lat * deg2rad;
  float cent_lonr__ = c_lon * deg2rad;
  float cent_latr__ = c_lat * deg2rad;

  float q1__ = cos(ob_lonr__ - cent_lonr__);
  float q2__ = cos(ob_latr__ - cent_latr__);
  float q3__ = cos(ob_latr__ + cent_latr__);

  float dij = (re * acos(0.5*((1.0+q1__)*q2__ - (1.0-q1__)*q3__) ) + 1.0);
  return dij;
}

// -----------------------------------------------------------------------------

void Gaussian_Thinning::applyFilter(const std::vector<bool> & apply,
                                    const oops::Variables & filtervars,
                                    std::vector<std::vector<bool>> & flagged) const {
  const size_t nobs = obsdb_.nlocs();
  const float re = Constants::mean_earth_rad;  // km
  const float deg2rad = Constants::deg2rad;
  const float half_circum = M_PI * re;
  const float grid_input_dist  = config_.getFloat("horizontal_mesh", 2 * M_PI * re / 360.0);
  const float vertical_mesh    = config_.getFloat("vertical_mesh", -99999.9);
                                                  // default in Pa; min value to thin over
  const float vertical_min     = config_.getFloat("vertical_min",    100.);
  const float vertical_max     = config_.getFloat("vertical_max", 110000.);

  ioda::ObsDataVector<float> odv_lat(obsdb_, "latitude", "MetaData");
  ioda::ObsDataVector<float> odv_lon(obsdb_, "longitude", "MetaData");
  ioda::ObsDataVector<float> odv_pres(obsdb_, "air_pressure", "MetaData");
  auto &lon = odv_lon[0];
  auto &lat = odv_lat[0];
  auto &pres = odv_pres[0];

  float grid_dist_km, grid_dist_deg, cur_grid_dist_deg;
  float cur_radius, cur_nlons;
  float clat, clon;
  float fn_vmesh;
  int nlon, i_mesh, i_vmesh;
  size_t n_lats, n_mesh, n_vmesh;
  float ob_lon, ob_lat, dist;
  int ob_idx, ob_vidx;

  oops::Log::debug() << "Gaussian_Thinning: config = " << config_ << std::endl;
  oops::Log::debug() << "Gaussian_Thinning: input thinning mesh = " << grid_input_dist << std::endl;

  // Balance input grid mesh size so that they are equally balanced around a sphere
  n_mesh = half_circum/grid_input_dist;
  grid_dist_km = half_circum/n_mesh;

  oops::Log::debug() << "Gaussian_Thinning: balanced thinning mesh = " << grid_dist_km << std::endl;

  // convert balanced mesh from km to degrees
  grid_dist_deg = grid_dist_km * 180.0 / half_circum;

  // get # of bins in lat direction
  n_lats = half_circum/grid_dist_km;
  std::vector<size_t> n_lons(n_lats);

  // construct bins in lon direction (variable as f(lat)); count total
  n_mesh = 0;
  clat = -90.0 + 0.5*grid_dist_deg;
  for (size_t i=0; i < n_lats; ++i) {   // iterate over n_lats
      /* n_lons is a function of n_lats; this calculates for each
         lat how many lon bins there are at the current lat.  
         The number of bins at lat scale to cos(clat) */
      n_lons[i] = cos(clat * deg2rad) * 360.0 / grid_dist_deg;
      // keep track of the total number of bins
      n_mesh += n_lons[i];
      // advance the latitude to the next bin
      clat += grid_dist_deg;
  }

  oops::Log::debug() << "Gaussian_Thinning: number of horizontal mesh bins = " << n_mesh
                                                                               << std::endl;

  // Create arrays of the bin centroid lat/lon
  std::vector<float> centroid_lat(nobs);
  std::vector<float> centroid_lon(nobs);

  i_mesh = 0;
  clat = -90.0 + 0.5*grid_dist_deg;
  for (size_t i=0; i < n_lats; ++i) {
     cur_grid_dist_deg = 360.0 / n_lons[i];
     clon = 0.5*cur_grid_dist_deg;
     for (size_t j=0; j < n_lons[i]; j++) {
        centroid_lat[i_mesh] = clat;
        centroid_lon[i_mesh] = clon;
        oops::Log::debug() << "Gaussian_Thinning: centroid_ll " << clon << ' ' << clat << std::endl;
        clon += cur_grid_dist_deg;
        i_mesh++;
     }
     clat += grid_dist_deg;
  }

  // calculate number of vertical bins
  if (vertical_mesh > 0.0) {
     fn_vmesh = (vertical_max - vertical_min) / vertical_mesh;
     n_vmesh = fn_vmesh;
     if ((fn_vmesh - n_vmesh) > 0.0) {++n_vmesh;}
  } else {
     n_vmesh = 1;
  }
  oops::Log::debug() << "Gaussian_Thinning: number of vertical mesh bins = " << n_vmesh
                                                                             << std::endl;

  /* winner and scoring arrays - each mesh point will have an index corresponding to the winner and the 
     score.  Those which win survive, those which do not will get the QCflags::thinned label slapped on
     them*/
  std::vector<std::vector<int>>   thin_idx(n_mesh, std::vector<int>(n_vmesh, -999));
  std::vector<std::vector<float>> score(n_mesh, std::vector<float>(n_vmesh, 99999.99));

  // loop through obs and test score
  for (size_t jobs = 0; jobs < nobs; ++jobs) {
     if (flags_[0][jobs] == 0) {     // only test obs that have passed QC to this point
        ob_lon = lon[jobs];
        ob_lat = lat[jobs];
        ob_idx = ll_to_idx(ob_lon, ob_lat, n_lats, n_lons);
        ob_vidx = pres_to_idx(pres[jobs], n_vmesh, vertical_mesh, vertical_max);

        dist = dist_to_centroid(ob_lon, ob_lat, centroid_lon[ob_idx], centroid_lat[ob_idx]);
        if (dist < score[ob_idx][ob_vidx]) {
           score[ob_idx][ob_vidx] = dist;
           thin_idx[ob_idx][ob_vidx] = jobs;
        }
        oops::Log::debug() << "Gaussian_Thinning: ob_ll_p " << lon[jobs] << " " << lat[jobs]
                                                            << " " << pres[jobs]<< std::endl;
     }
  }

  std::vector<bool> thin(nobs, true);

  // loop through scoring mesh and project obs that pass thinning from mesh to obs space
  for (size_t i=0; i < n_mesh; i++) {
     for (size_t k=0; k < n_vmesh; k++) {
        if (thin_idx[i][k] > 0) {
           oops::Log::debug() << "Gaussian_Thinning: selob_ll " << lon[thin_idx[i][k]] << " "
                                                                << lat[thin_idx[i][k]] << " "
                                                                << pres[thin_idx[i][k]]
                                                                << std::endl;
           thin[thin_idx[i][k]] = false;
        }
     }
  }

  // project the QC across all varialbes and fail all obs that do not pass thinning
  for (size_t jv = 0; jv < filtervars.size(); ++jv) {
     for (size_t jobs = 0; jobs < nobs; ++jobs) {
        if ( apply[jobs] && thin[jobs] ) flagged[jv][jobs] = true;
     }
  }
}

// -----------------------------------------------------------------------------

void Gaussian_Thinning::print(std::ostream & os) const {
  os << "Gaussian_Thinning: config = " << config_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
