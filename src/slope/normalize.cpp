#include "normalize.h"

namespace slope {

JitNormalization
normalize(Eigen::MatrixXd& x,
          Eigen::VectorXd& x_centers,
          Eigen::VectorXd& x_scales,
          const std::string& centering_type,
          const std::string& scaling_type,
          const bool modify_x)
{
  const int p = x.cols();

  computeCenters(x_centers, x, centering_type);
  computeScales(x_scales, x, scaling_type);

  if ((x_scales.array().abs() == 0.0).any()) {
    throw std::invalid_argument("One or more columns have zero variance");
  }

  bool center = centering_type != "none";
  bool scale = scaling_type != "none";
  bool center_jit = center && !modify_x;
  bool scale_jit = scale && !modify_x;

  JitNormalization jit_normalization;

  if (center_jit && scale_jit) {
    jit_normalization = JitNormalization::Both;
  } else if (center_jit) {
    jit_normalization = JitNormalization::Center;
  } else if (scale_jit) {
    jit_normalization = JitNormalization::Scale;
  } else {
    jit_normalization = JitNormalization::None;
  }

  if (modify_x && (center || scale)) {
    for (int j = 0; j < p; ++j) {
      if (center) {
        x.col(j).array() -= x_centers(j);
      }
      if (scale) {
        x.col(j).array() /= x_scales(j);
      }
    }
  }

  return jit_normalization;
}

JitNormalization
normalize(Eigen::SparseMatrix<double>& x,
          Eigen::VectorXd& x_centers,
          Eigen::VectorXd& x_scales,
          const std::string& centering_type,
          const std::string& scaling_type,
          const bool)
{
  computeCenters(x_centers, x, centering_type);
  computeScales(x_scales, x, scaling_type);

  bool center = centering_type != "none";
  bool scale = scaling_type != "none";
  bool center_jit = center;
  bool scale_jit = scale;

  JitNormalization jit_normalization;

  if (center_jit && scale_jit) {
    jit_normalization = JitNormalization::Both;
  } else if (center_jit) {
    jit_normalization = JitNormalization::Center;
  } else if (scale_jit) {
    jit_normalization = JitNormalization::Scale;
  } else {
    jit_normalization = JitNormalization::None;
  }

  // TODO: Implement in-place scaling for sparse matrices.
  // if (modify_x && scaling_type != "none") {
  //   for (int j = 0; j < x.cols(); ++j) {
  //     for (Eigen::SparseMatrix<double>::InnerIterator it(x, j); it; ++it) {
  //       it.valueRef() = it.value() / x_scales(j);
  //     }
  //   }
  // }

  return jit_normalization;
}

// TODO: Make this work for sparse beta
std::tuple<Eigen::VectorXd, Eigen::MatrixXd>
rescaleCoefficients(const Eigen::VectorXd& beta0,
                    const Eigen::SparseMatrix<double>& beta,
                    const Eigen::VectorXd& x_centers,
                    const Eigen::VectorXd& x_scales,
                    const bool intercept)
{
  int m = beta0.size();

  bool centering = x_centers.size() > 0;
  bool scaling = x_scales.size() > 0;

  Eigen::MatrixXd beta0_out = beta0;
  Eigen::SparseMatrix<double> beta_out = beta;

  if (centering || scaling) {
    for (int k = 0; k < m; ++k) {
      double x_bar_beta_sum = 0.0;

      for (Eigen::SparseMatrix<double>::InnerIterator it(beta_out, k); it;
           ++it) {
        int j = it.row();

        if (scaling) {
          it.valueRef() /= x_scales(j);
          // beta_out(j, k) /= x_scales(j);
        }
        if (centering) {
          x_bar_beta_sum += x_centers(j) * it.valueRef();
        }
      }

      if (intercept) {
        beta0_out(k) -= x_bar_beta_sum;
      }
    }
  }

  return { beta0_out, beta_out };
}

} // namespace slope
