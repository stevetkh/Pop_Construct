#include <TMB.hpp>                                

template<class Type>
Type objective_function<Type>::operator() ()
{
  using Eigen::Matrix;
  using Eigen::Map;
  using Eigen::Dynamic;
  using namespace density;

  typedef Matrix<Type, Dynamic, Dynamic> MatrixXXt;
  typedef Map<Matrix<Type, Dynamic, 1>> MapVectorXt;
  typedef Map<MatrixXXt> MapMatrixXXt;

  DATA_MATRIX(LQ_baseline_f);
  PARAMETER_VECTOR(log_D_f);
  PARAMETER_VECTOR(log_E_f);
  PARAMETER_VECTOR(log_F_f);

  PARAMETER_VECTOR(h_params_f);
  PARAMETER_VECTOR(k_params_f);

  DATA_VECTOR(DEF_log_age);

  vector<Type> log_mort_rates_0_f = h_params_f - log(1.0 - 0.5 * exp(h_params_f)) + log(1.0 / 5.0);
  vector<Type> h2_params_f = h_params_f * h_params_f;

  matrix<Type> LQ_params_mat_f(4, h_params_f.size());
  for(int i=0; i < h_params_f.size(); i++){
  LQ_params_mat_f(0,i) = Type(1.0);
  }
  LQ_params_mat_f.row(1) = h_params_f;
  LQ_params_mat_f.row(2) = h2_params_f;
  LQ_params_mat_f.row(3) = k_params_f;

  matrix<Type> log_mx_mat_f(LQ_baseline_f.rows() + 1, h_params_f.size());
  log_mx_mat_f.row(0) = log_mort_rates_0_f;
  log_mx_mat_f.bottomRows(LQ_baseline_f.rows()) = LQ_baseline_f*LQ_params_mat_f;

  matrix<Type> mx_mat_f = exp(log_mx_mat_f.array());
  matrix<Type> log_DEF_mat_f;
  for(int i=0; i < h_params_f.size(); i++){
    log_DEF_mat_f.col(i) = log_D_f(i) - exp(log_E_f(i)) * (DEF_log_age - log_F_f(i)) * (DEF_log_age - log_F_f(i));
  }
  matrix<Type> DEF_mat_f = exp(log_DEF_mat_f.array());
}


