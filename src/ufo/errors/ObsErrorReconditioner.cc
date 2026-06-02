/*
 * (C) Crown copyright 2025 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#include <vector>

#include "ufo/errors/ObsErrorReconditioner.h"

namespace ufo {

constexpr char ObsErrorReconditionerMethodParameterTraitsHelper::enumTypeName[];
constexpr util::NamedEnumerator<ObsErrorReconditionerMethod>
ObsErrorReconditionerMethodParameterTraitsHelper::namedValues[];

ObsErrorReconditioner::ObsErrorReconditioner(const Parameters_ & params)
  : params_(params) {
  // Checking valid reconditioning options if reconditioning specified.
  ufo::ObsErrorReconditionerMethod recon_method = params_.ReconMethod.value();
  size_t nvalid_options = 0;
  switch (recon_method) {
    case ufo::ObsErrorReconditionerMethod::MINIMUMEIGENVALUE:
      nvalid_options += static_cast<int>(params_.kFrac.value() != boost::none);
      nvalid_options += static_cast<int>(params_.Threshold.value() != boost::none);
      if (nvalid_options == 0) {
        throw eckit::BadParameter("No viable reconditioning metric"
                                  " for minimum eigenvalue provided.",
                                  Here());
      } else if (nvalid_options > 1) {
        throw eckit::BadParameter("Too many reconditioning metrics provided.", Here());
      }
      break;
    case ufo::ObsErrorReconditionerMethod::RIDGEREGRESSION:
      nvalid_options += static_cast<int>(params_.kFrac.value() != boost::none);
      nvalid_options += static_cast<int>(params_.Shift.value() != boost::none);
      if (nvalid_options == 0) {
        throw eckit::BadParameter("No viable reconditioning metric"
                                  " for ridge regression provided.",
                                  Here());
      } else if (nvalid_options > 1) {
        throw eckit::BadParameter("Too many reconditioning metrics provided.", Here());
      }
      break;
    case ufo::ObsErrorReconditionerMethod::NORECONDITIONING:
      oops::Log::trace() << "'No reconditioning' option selected, "
                            "recondition method can be tested, "
                            "R matrix should not change.\n";
      break;
  }
}  // ObsErrorReconditioner::ObsErrorReconditioner

void ObsErrorReconditioner::recondition(Eigen::MatrixXd & R) const {
    // Get reconditioner method and if no reconditioning then exit straight away
    // this avoids unneccessary computation.
    ufo::ObsErrorReconditionerMethod recon_method = params_.ReconMethod.value();
    if (recon_method == ufo::ObsErrorReconditionerMethod::NORECONDITIONING) {
      oops::Log::trace() << "No reconditioning selected, skipping reconditioning.\n";
      return;
    }

    // Check square matrix
    size_t nrows = R.rows();
    size_t ncols = R.cols();
    assert(nrows == ncols);

    // Performing eigendecomposition
    oops::Log::trace() << "R before reconditioning:\n" << R << std::endl << std::endl;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;
    solver.compute(R);
    Eigen::VectorXd evals = solver.eigenvalues();
    const Eigen::MatrixXd & evecs = solver.eigenvectors();

    // Ensure positive definiteness
    double eval_min = evals.minCoeff();
    if (eval_min <= 0.0) {
      oops::Log::debug() << "eval_min = " << eval_min
                         << " , performing a ridge regression to ensure positive definiteness\n";
      double alpha = 1.0 + 1e-15;
      alpha = (eval_min == 0.0) ? alpha - 1.0 : alpha * abs(eval_min);
      for (size_t jvar = 0; jvar < nrows; ++jvar) {
          evals[jvar] += alpha;
      }
    }
    const double precondno = evals.maxCoeff() / evals.minCoeff();
    oops::Log::trace() << "Condition no. before = " << precondno << std::endl;

    // Reconditioning
    double threshold = 0.0;
    double delta = 0.0;
    switch (recon_method) {
      case ufo::ObsErrorReconditionerMethod::MINIMUMEIGENVALUE:
        oops::Log::trace() << "Performing Minimum Eigen Value reconditioning\n";
        // Determining threshold from options
        if (params_.Threshold.value() != boost::none) {
          threshold = params_.Threshold.value().value();
        } else if (params_.kFrac.value() != boost::none) {
          const double kmax = precondno * params_.kFrac.value().value();
          if (kmax < 1.0) {
            oops::Log::warning() << "Unviable kmax = "
                                 << kmax
                                 << ", skipping reconditioning\n";
            return;
          }
          threshold = evals.maxCoeff() / kmax;
        }
        // Adjusting eigenvalues
        for (size_t jvar = 0; jvar < nrows; jvar++) {
          if (evals[jvar] <= threshold) {
            evals[jvar] = threshold;
          }
        }
        break;
      case ufo::ObsErrorReconditionerMethod::RIDGEREGRESSION:
        oops::Log::trace() << "Performing Ridge Regression reconditioning\n";
        // Determining delta from options
        if (params_.kFrac.value() != boost::none) {
          const double kmax = precondno * params_.kFrac.value().value();
          if (kmax < 1.0) {
            oops::Log::warning() << "Unviable kmax = "
                                 << kmax
                                 << ", skipping reconditioning\n";
            return;
          }
          delta = (evals.maxCoeff() - evals.minCoeff() * kmax)/(kmax - 1);
        } else if (params_.Shift.value() != boost::none) {
          delta = params_.Shift.value().value();
        }
        // Adjusting eigenvalues
        for (size_t jvar = 0; jvar < nrows; ++jvar) {
          evals[jvar] += delta;
        }
        break;
      case ufo::ObsErrorReconditionerMethod::NORECONDITIONING:
        break;
    }

    // Fail-safe for if reconditioning
    // produces non positive definite matrix
    eval_min = evals.minCoeff();
    if (eval_min <= 0.0) {
      oops::Log::trace() << "eval_min = " << eval_min
                         << " , skipping reconditioning\n";
      return;
    }

    // Re-evaluating R
    R = evecs * evals.asDiagonal() * evecs.transpose();
    oops::Log::trace() << "R after reconditioning:\n" << R << std::endl << std::endl;
    const double condno = evals.maxCoeff() / evals.minCoeff();
    oops::Log::trace() << "Condition no. after = " << condno << std::endl;
    oops::Log::trace() << "Ratio of condition numbers = " << condno / precondno << std::endl;
}  // ObsErrorReconditioner::recondition

}  // namespace ufo
