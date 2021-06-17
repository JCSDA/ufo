/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <Eigen/Dense>
#include <algorithm>
#include <set>

#include "ioda/ObsDataVector.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/missingValues.h"
#include "ufo/filters/obsfunctions/CloudCostFunction.h"
#include "ufo/filters/Variable.h"
#include "ufo/utils/ufo_utils.interface.h"

namespace ufo {

static ObsFunctionMaker<CloudCostFunction> makerCloudCostFunction_("CloudCostFunction");

CloudCostFunction::CloudCostFunction(const eckit::LocalConfiguration & conf)
  : invars_() {
  // Initialize options
  options_.deserialize(conf);

  // List of field names for B matrix
  fields_ = options_.field_names.value();

  // Get channels for computing scattering index from options
  std::set<int> chanset = oops::parseIntSet(options_.chanlist.value());
  channels_.assign(chanset.begin(), chanset.end());

  // List of required data
  for (size_t i = 0; i < fields_.size(); ++i) {
    invars_ += Variable("brightness_temperature_jacobian_"+fields_[i]+"@ObsDiag", channels_);
  }
  invars_ += Variable("brightness_temperature@ObsValue", channels_);
  invars_ += Variable("brightness_temperature@"+options_.HofXGroup.value(), channels_);
  invars_ += Variable("latitude@MetaData");
  if (options_.qtotal_lnq_gkg.value()) {
    invars_ += Variable("specific_humidity@GeoVaLs");
    invars_ += Variable("mass_content_of_cloud_liquid_water_in_atmosphere_layer@GeoVaLs");
    invars_ += Variable("mass_content_of_cloud_ice_in_atmosphere_layer@GeoVaLs");
    invars_ += Variable("air_pressure@GeoVaLs");
    invars_ += Variable("air_temperature@GeoVaLs");
    invars_ += Variable("surface_pressure@GeoVaLs");
    invars_ += Variable("surface_temperature@GeoVaLs");
    invars_ += Variable("specific_humidity_at_two_meters_above_surface@GeoVaLs");
  }
}

// -----------------------------------------------------------------------------

void CloudCostFunction::compute(const ObsFilterData & in,
                                    ioda::ObsDataVector<float> & out) const {
  // Get dimensions
  const size_t nlocs = in.nlocs();
  const size_t nchans = channels_.size();
  ASSERT(nchans > 0);
  ASSERT(out.nvars() == 1);

  // B, R error covariance objects
  eckit::LocalConfiguration bMatrixConf;
  bMatrixConf.set("BMatrix", options_.bmatrix_filepath.value());
  bMatrixConf.set("background fields", options_.field_names.value());
  bMatrixConf.set("qtotal", options_.qtotal_lnq_gkg.value());
  MetOfficeBMatrixStatic staticB(bMatrixConf);

  eckit::LocalConfiguration rMatrixConf;
  rMatrixConf.set("RMatrix", options_.rmatrix_filepath.value());
  MetOfficeRMatrixRadiance staticR(rMatrixConf);

  bool split_rain = options_.qtotal_split_rain.value();

  const std::string clw_name = "mass_content_of_cloud_liquid_water_in_atmosphere_layer";
  const std::string ciw_name = "mass_content_of_cloud_ice_in_atmosphere_layer";
  std::vector<float> gv_pres(nlocs), gv_temp(nlocs), gv_qgas(nlocs), gv_clw(nlocs), gv_ciw(nlocs),
                     humidity_total(nlocs);

  // Determine if pressure is ascending or descending (B-matrix assumption)
  size_t np = in.nlevs(Variable("air_pressure@GeoVaLs"));
  std::vector<float> gv_pres_1(nlocs), gv_pres_N(nlocs);
  in.get(Variable("air_pressure@GeoVaLs"), 0, gv_pres_1);
  in.get(Variable("air_pressure@GeoVaLs"), np - 1, gv_pres_N);
  const float missing = util::missingValue(missing);
  ASSERT(gv_pres_1[0] != missing);
  ASSERT(gv_pres_N[0] != missing);
  bool p_ascending = (gv_pres_N[0] > gv_pres_1[0]);

  // Assemble combined Jacobian from component fields
  std::vector<std::vector<std::vector<float>>>
       jac_vec(nlocs, std::vector<std::vector<float>>(nchans, std::vector<float>()));
  for (size_t ifield = 0; ifield < fields_.size(); ++ifield) {
    if (options_.qtotal_lnq_gkg.value() &&
            (fields_[ifield] == clw_name || fields_[ifield] == ciw_name)) {
      // qtotal ln(g/kg) Jacobian calculated when "specific_humidity" is reached in field list
      continue;
    }
    std::string jac_name = "brightness_temperature_jacobian_"+fields_[ifield]+"@ObsDiag";
    size_t nlevs = in.nlevs(Variable(jac_name, channels_)[0]);
    std::vector<float> jac_store(nlocs);
    for (size_t ilev = 0; ilev < nlevs; ++ilev) {
      const int level_gv = (p_ascending ? ilev : nlevs-ilev-1);
      const int level_jac = (options_.reverse_Jacobian.value() ? nlevs-level_gv-1 : level_gv);
      if (fields_[ifield] == "specific_humidity" && options_.qtotal_lnq_gkg.value()) {
        in.get(Variable("air_pressure@GeoVaLs"), level_gv, gv_pres);
        in.get(Variable("air_temperature@GeoVaLs"), level_gv, gv_temp);
        in.get(Variable("specific_humidity@GeoVaLs"), level_gv, gv_qgas);
        in.get(Variable(clw_name+"@GeoVaLs"), level_gv, gv_clw);
        in.get(Variable(ciw_name+"@GeoVaLs"), level_gv, gv_ciw);
        std::vector<float> qsaturated(nlocs);
        ufo_ops_satrad_qsatwat_f90(qsaturated.data(), gv_temp.data(), gv_pres.data(),
                                   static_cast<int>(nlocs));
        for (size_t iloc = 0; iloc < nlocs; ++iloc) {
          // Ensure specific humidity is within limits
          gv_qgas[iloc] = std::max(gv_qgas[iloc], options_.min_q.value());
          gv_qgas[iloc] = std::min(gv_qgas[iloc], qsaturated[iloc]);
          humidity_total[iloc] = gv_qgas[iloc] + gv_clw[iloc] + gv_ciw[iloc];
        }
      }

      for (size_t ichan = 0; ichan < nchans; ++ichan) {
        in.get(Variable(jac_name, channels_)[ichan], level_jac, jac_store);

        if (fields_[ifield] == "specific_humidity" && options_.qtotal_lnq_gkg.value()) {
          std::vector<float> jac_clw(nlocs), jac_ciw(nlocs);
          in.get(Variable("brightness_temperature_jacobian_"+clw_name+"@ObsDiag", channels_)[ichan],
                 level_jac, jac_clw);
          in.get(Variable("brightness_temperature_jacobian_"+ciw_name+"@ObsDiag", channels_)[ichan],
                 level_jac, jac_ciw);
          std::vector<float> dq_dqtotal(nlocs), dql_dqtotal(nlocs), dqi_dqtotal(nlocs);
          int qsplit_mode = 2;  // compute derivatives
          ufo_ops_satrad_qsplit_f90(qsplit_mode, static_cast<int>(nlocs), gv_pres.data(),
                                    gv_temp.data(), humidity_total.data(), dq_dqtotal.data(),
                                    dql_dqtotal.data(), dqi_dqtotal.data(), split_rain);
          // Jacobian dy/dx for observation y, humdity x in units kg/kg
          // For alternative units of B-matrix humidity z in ln(g/kg)
          // chain rule gives dy/dz = x.(dy/dx)
          // Gradient due to ice is ignored unless we are using scattering radiative transfer
          // dTb/dln(qt) = qt*(dTb/dq*dq/dqt + dTb/dql*dql/dqt [+ dTb/dqi*dqi/dqt])
          for (size_t iloc = 0; iloc < nlocs; ++iloc) {
            jac_store[iloc] *= dq_dqtotal[iloc];
            jac_store[iloc] += jac_clw[iloc]*dql_dqtotal[iloc];
            if (options_.scattering_switch.value()) {
              jac_store[iloc] += jac_ciw[iloc]*dqi_dqtotal[iloc];
            }
            jac_store[iloc] *= humidity_total[iloc];
          }
        }

        if (fields_[ifield] == "specific_humidity_at_two_meters_above_surface"
            && options_.qtotal_lnq_gkg.value()) {
          in.get(Variable("surface_pressure@GeoVaLs"), level_gv, gv_pres);
          in.get(Variable("surface_temperature@GeoVaLs"), level_gv, gv_temp);
          in.get(Variable("specific_humidity_at_two_meters_above_surface@GeoVaLs"),
                 level_gv, gv_qgas);
          std::vector<float> qsaturated(nlocs);
          ufo_ops_satrad_qsatwat_f90(qsaturated.data(), gv_temp.data(), gv_pres.data(),
                                     static_cast<int>(nlocs));
          for (size_t iloc = 0; iloc < nlocs; ++iloc) {
            gv_qgas[iloc] = std::max(gv_qgas[iloc], options_.min_q.value());
            gv_qgas[iloc] = std::min(gv_qgas[iloc], qsaturated[iloc]);
            // dTb/dln(q2m) = q2m*(dTb/dq2m)
            jac_store[iloc] *= gv_qgas[iloc];
          }
        }

        for (size_t iloc = 0; iloc < nlocs; ++iloc) {
          jac_vec[iloc][ichan].push_back(jac_store[iloc]);
        }
      }
    }
  }

  const size_t sizeB = staticB.getsize();
  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    for (size_t ichan = 0; ichan < nchans; ++ichan) {
      ASSERT(jac_vec[iloc][ichan].size() == sizeB);
    }
  }

  // Get departures = ObsValue - HofX (where HofX is bias corrected)
  std::vector<std::vector<float>> departures(nlocs, std::vector<float>(nchans));
  std::vector<float> obsvalues(nlocs);
  std::vector<float> bgvalues(nlocs);
  std::vector<bool> is_out_of_bounds(nlocs, false);
  for (size_t ichan = 0; ichan < nchans; ++ichan) {
    in.get(Variable("brightness_temperature@ObsValue", channels_)[ichan], obsvalues);
    in.get(Variable("brightness_temperature@"+options_.HofXGroup.value(), channels_)[ichan],
           bgvalues);
    for (size_t iloc = 0; iloc < nlocs; ++iloc) {
      departures[iloc][ichan] = obsvalues[iloc] - bgvalues[iloc];
      // Flag observations outside expected bounds
      if (obsvalues[iloc] < options_.minTb.value() || obsvalues[iloc] > options_.maxTb.value()) {
        is_out_of_bounds[iloc] = true;
      }
    }
  }

  std::vector<float> latitude(nlocs);
  in.get(Variable("latitude@MetaData"), latitude);
  Eigen::MatrixXf Hmatrix(nchans, sizeB);

  for (size_t iloc = 0; iloc < nlocs; ++iloc) {
    if (is_out_of_bounds[iloc]) {
      out[0][iloc] = options_.maxCost.value();
      continue;
    }

    for (size_t ichan = 0; ichan < nchans; ++ichan) {
      Hmatrix.row(ichan) = Eigen::Map<Eigen::VectorXf>(jac_vec[iloc][ichan].data(), sizeB);
    }

    // Matrix of departures dy
    Eigen::Map<Eigen::VectorXf> dy(departures[iloc].data(), nchans);

    // Calculate Scratch_matrix = H.B.H^T + R
    Eigen::MatrixXf BHT;
    staticB.multiply(latitude[iloc], Hmatrix.transpose(), BHT);
    Eigen::MatrixXf HBHT = Hmatrix*BHT;
    Eigen::MatrixXf Scratch_matrix;
    staticR.add(channels_, HBHT, Scratch_matrix);

    // Calculate Scratch_matrix2 = Scratch_matrix^-1.dy using Cholesky decomposition
    Eigen::LLT<Eigen::MatrixXf> decomposition(Scratch_matrix);
    if (decomposition.info() == Eigen::NumericalIssue) {
      oops::Log::warning() <<
        "CloudCostFunction Scratch_matrix appears not to be positive definite" << std::endl;
      out[0][iloc] = options_.maxCost.value();
      continue;
    }
    Eigen::VectorXf Scratch_matrix2 = decomposition.solve(dy);

    // Final cost
    float Cost_final = 0.5*dy.transpose()*Scratch_matrix2;
    Cost_final /= static_cast<float>(nchans);  // normalise by number of channels
    out[0][iloc] = std::min(Cost_final, options_.maxCost.value());
  }
}

// -----------------------------------------------------------------------------

const ufo::Variables & CloudCostFunction::requiredVariables() const {
  return invars_;
}

// -----------------------------------------------------------------------------

}  // namespace ufo
