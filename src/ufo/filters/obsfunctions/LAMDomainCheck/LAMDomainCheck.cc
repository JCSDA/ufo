/*
 * (C) Copyright 2020 NOAA NWS NCEP EMC
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "ufo/filters/obsfunctions/LAMDomainCheck/LAMDomainCheck.h"

#include "eckit/exception/Exceptions.h"

#include "oops/util/missingValues.h"

#include "ufo/filters/Variable.h"

namespace ufo {

static ObsFunctionMaker<LAMDomainCheck> makerObsFuncLAMDomainCheck_("LAMDomainCheck");

// -----------------------------------------------------------------------------

LAMDomainCheck::LAMDomainCheck(const eckit::LocalConfiguration & conf)
  : invars_() {
  oops::Log::debug() << "LAMDomainCheck: config = " << conf << std::endl;
  // Initialize options
  options_.deserialize(conf);

  // We must know the latitude of each observation
  invars_ += Variable("MetaData/latitude");
  // We must know the longitude of each observation
  invars_ += Variable("MetaData/longitude");
}

// -----------------------------------------------------------------------------

LAMDomainCheck::~LAMDomainCheck() {}

// -----------------------------------------------------------------------------
/*! \brief LAMDomainCheck::compute
*
* \details The LAMDomainCheck is an obsfunction to compute if an observation is
* located inside a specified limited area model domain.
* The LAMDomainCheck obsfunction returns a value of 1 if the observation is determined
* to be located inside the specified domain and a value of 0 if it lies outside.
*
* In the UFO obs functions YAML, first one must define the map_projection.
* The following values of map_projection are currently supported:
* * "gnomonic_ed" - the ESG grid used by FV3-LAM
*
* The option 'save: true' will save the computed value to the output IODA file as
* 'DerivedValue/LAMDomainCheck' (default is false).
*
* The additional parameters to be defined in the options section
* of the obs function YAML will depend on the choice of map_projection used.
* For gnomonic_ed:
* * a - ESG alpha parameter (default 0)
* * k - ESG kappa parameter (default 0)
* * plat - ESG center point latitude (degrees; default 0)
* * plon - ESG center point longitude (degrees; default 0)
* * pazi - ESG azimuthal angle (radians; default 0)
* * dx - grid spacing in x (degrees; default 1)
* * dy - grid spacing in y (degrees, default 1)
* * npx - number of gridpoints in x (default 2)
* * npy - number of gridpoints in y (default 2)
* * nbdy - number of gridpoints from lateral boundary where observations are rejected in this buffer zone (default 0)
* If using FV3-LAM, the above values are provided in the netCDF grid file attributes,
* but this option will work for any regional model that utilizes the
* Extended Schmidt Gnomonic grid developed by R. Jim Purser:
* https://dtcenter.org/sites/default/files/events/2020/2-purser-james.pdf
* "The Extended Schmidt Gnomonic grid for regional applications"
* by R. J. Purser, D. Jovic, G. Ketefian, T. Black, J. Beck, J. Dong, J. Carley.
* UFS Users' Workshop, July 27--29, 2020.
*
*/

void LAMDomainCheck::compute(const ObsFilterData & in,
                                  ioda::ObsDataVector<float> & out) const {
  const size_t nlocs = in.nlocs();
  const float missing = util::missingValue<float>();

  // Ensure that only one output variable is expected.
  ASSERT(out.nvars() == 1);

  std::vector<int> iidx(nlocs);

  // Retrieve the latitude and longitude.
  std::vector<float> latitude;
  std::vector<float> longitude;
  in.get(Variable("MetaData/latitude"), latitude);
  in.get(Variable("MetaData/longitude"), longitude);

  // get options based off the name of the map projection
  if (options_.mapproj.value() == "gnomonic_ed") {
    // ESG used in FV3-LAM
    const float a = options_.esg_a.value();
    const float k = options_.esg_k.value();
    const float plat = options_.esg_plat.value();
    const float plon = options_.esg_plon.value();
    const float pazi = options_.esg_pazi.value();
    const float dx = options_.esg_dx.value();
    const float dy = options_.esg_dy.value();
    const int npx = options_.esg_npx.value();
    const int npy = options_.esg_npy.value();
    const int nbdy = options_.esg_nbdy.value();
    for (size_t jj = 0; jj < nlocs; ++jj) {
      lam_domaincheck_esg_f90(a, k, plat, plon, pazi, npx, npy,
                         dx, dy, nbdy, latitude[jj], longitude[jj], iidx[jj]);
      out[0][jj] = static_cast<float>(iidx[jj]);
    }
  } else if (options_.mapproj.value() == "circle") {
    const float cenlat = options_.cenlat.value();
    const float cenlon = options_.cenlon.value();
    const float radius = options_.radius.value();
    for (size_t jj = 0; jj < nlocs; ++jj) {
      // calculate great-circle distance on sphere
      lam_domaincheck_circle_f90(cenlat, cenlon, radius,
                         latitude[jj], longitude[jj], iidx[jj]);
      out[0][jj] = static_cast<float>(iidx[jj]);
    }
  } else {
    // throw exception for unsupported projection
    std::string errString = " is not a supported map projection. Fatal error!!!";
    oops::Log::error() << options_.mapproj.value() << errString;
    throw eckit::BadValue(errString);
  }

  if (options_.save) {
    out.save("DerivedValue");
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & LAMDomainCheck::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
