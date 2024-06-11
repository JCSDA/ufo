/*
 * (C) Copyright 2020-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "ufo/operators/marine/chleuzintegr/ObsChlEuzIntegr.h"

#include <ostream>
#include <string>
#include <vector>

#include "ioda/ObsVector.h"

#include "oops/base/Variables.h"

#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/utils/Constants.h"

namespace ufo {

// -----------------------------------------------------------------------------
static ObsOperatorMaker<ObsChlEuzIntegr> makerChlEuzIntegr_("Chlorophyll Ocean Color");
// -----------------------------------------------------------------------------

ObsChlEuzIntegr::ObsChlEuzIntegr(const ioda::ObsSpace & odb,
                                 const Parameters_ & params)
  : ObsOperatorBase(odb), varin_()
{
  const std::vector<std::string> vvin{"mass_concentration_of_chlorophyll_in_sea_water",
                                      "sea_water_cell_thickness"};
  varin_.reset(new oops::Variables(vvin));
  oops::Log::trace() << "ObsChlEuzIntegr created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsChlEuzIntegr::~ObsChlEuzIntegr() {
  oops::Log::trace() << "ObsChlEuzIntegr destructed" << std::endl;
}

// -----------------------------------------------------------------------------
void ObsChlEuzIntegr::simulateObs(const GeoVaLs & gv, ioda::ObsVector & ovec,
                                  ObsDiagnostics & d, const QCFlags_t & qc_flags) const {
  int nlocs = ovec.size();
  int nlevs = gv.nlevs(oops::Variable{"mass_concentration_of_chlorophyll_in_sea_water"});

  const double missing = util::missingValue<double>();

  // common vectors storage
  std::vector <double> tmp(nlocs, 0.0);

  // Retrieve the chlorophyll and cell thickness
  std::vector<std::vector<double>> chl;
  std::vector<std::vector<double>> h;
  for ( std::size_t k = 0; k < nlevs; ++k ) {
    gv.getAtLevel(tmp, oops::Variable{"sea_water_cell_thickness"}, k);
    h.push_back(tmp);
    gv.getAtLevel(tmp, oops::Variable{"mass_concentration_of_chlorophyll_in_sea_water"}, k);
    chl.push_back(tmp);
  }

  // Calculate mean chlorophyll averaged over euphotic layer (euz_mod)
  for ( std::size_t i = 0; i < nlocs; ++i ) {
    // skip missing geovals
    if (chl[0][i] == missing) {
      ovec[i] = missing;
      continue;
    }

    double euz = Constants::euzc_0 * pow(chl[0][i], Constants::euzc_1);
    double euz_mod = 0.0;
    int elev = 0;
    for ( std::size_t k = 0; k < nlevs; ++k ) {
      if (euz_mod < euz) {
        euz_mod += h[k][i];
        elev++;
      }
    }
    ovec[i] = 0.0;
    for ( std::size_t k = 0; k < elev; ++k ) {
      ovec[i] += chl[k][i] * h[k][i] / euz_mod;
    }
  }
  oops::Log::trace() << "ObsChlEuzIntegr: observation operator run" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsChlEuzIntegr::print(std::ostream & os) const {
  os << "Chlorophyll Ocean Color obs operator";
}

// -----------------------------------------------------------------------------

}  // namespace ufo
