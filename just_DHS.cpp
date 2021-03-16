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
 
  DATA_VECTOR(LQ_baseline_mx_DX_f);
  DATA_SPARSE_MATRIX(h_DX_f);
  DATA_SPARSE_MATRIX(h2_DX_f);
  DATA_SPARSE_MATRIX(k_DX_f);
  DATA_SPARSE_MATRIX(tp_DX_f);
  DATA_SPARSE_MATRIX(penal_tp);
  DATA_SPARSE_MATRIX(penal_tp_0);
  DATA_SPARSE_MATRIX(null_penal_tp);
  DATA_MATRIX(LQ_baseline_f);
  DATA_VECTOR(df);
  DATA_VECTOR(Ef);

  DATA_VECTOR(LQ_baseline_mx_DX_m);
  DATA_SPARSE_MATRIX(h_DX_m);
  DATA_SPARSE_MATRIX(h2_DX_m);
  DATA_SPARSE_MATRIX(k_DX_m);
  DATA_SPARSE_MATRIX(tp_DX_m);
  DATA_MATRIX(LQ_baseline_m);
  DATA_VECTOR(dm);
  DATA_VECTOR(Em);

  PARAMETER_VECTOR(log_marginal_var_h);
  PARAMETER_VECTOR(log_marginal_var_k);
  PARAMETER(log_lambda_tp);
  PARAMETER(log_lambda_tp_0_inflated_sd);
  PARAMETER_VECTOR(log_dispersion);
 


  PARAMETER_VECTOR(h_params_f);
  PARAMETER_VECTOR(k_params_f);
  PARAMETER_VECTOR(logit_rho_h);
  PARAMETER_VECTOR(logit_rho_k);
  PARAMETER(k_constant_f);
  PARAMETER_VECTOR(tp_params);
  vector<Type> h2_params_f = h_params_f*h_params_f;


  PARAMETER_VECTOR(h_params_m);
  PARAMETER_VECTOR(k_params_m);
  PARAMETER(k_constant_m);
  vector<Type> h2_params_m = h_params_m*h_params_m;

  Type nll(0.0);

  nll -= dnorm(log_marginal_var_h, Type(0.0), Type(5.0), 1).sum();
  nll -= dnorm(log_marginal_var_k, Type(0.0), Type(5.0), 1).sum();

  nll -= dnorm(logit_rho_h, Type(0.0), Type(1.0), 1).sum();
  nll -= dnorm(logit_rho_k, Type(0.0), Type(1.0), 1).sum();

  nll -= dnorm(k_constant_f, Type(0.0), Type(1.0), 1);
  nll -= dnorm(k_constant_m, Type(0.0), Type(1.0), 1);

  nll -= dnorm(log_lambda_tp, Type(0.0), Type(5.0), 1);
  nll -= dnorm(log_lambda_tp_0_inflated_sd, Type(0.0), Type(5.0), 1);
  
  nll -= dnorm(log_dispersion, Type(0.0), Type(5.0), 1).sum();

  vector<Type> h_innov_f, h_innov_m;
  h_innov_f = diff(h_params_f);
  h_innov_m = diff(h_params_m);

  nll += SCALE(AR1(invlogit(logit_rho_h(0))), exp(0.5*log_marginal_var_h(0)))(h_innov_f);
  nll += SCALE(AR1(invlogit(logit_rho_h(1))), exp(0.5*log_marginal_var_h(1)))(h_innov_m);
  nll -= dnorm(h_params_f(0), Type(-2.0), Type(5.0), 1);
  nll -= dnorm(h_params_m(0), Type(-2.0), Type(5.0), 1);

  nll += SCALE(AR1(invlogit(logit_rho_k(0))), exp(0.5*log_marginal_var_k(0)))(k_params_f - k_constant_f);
  nll += SCALE(AR1(invlogit(logit_rho_k(1))), exp(0.5*log_marginal_var_k(1)))(k_params_m - k_constant_m);

  SparseMatrix<Type> QQ_tp = exp(log_lambda_tp)*penal_tp + exp(-2*log_lambda_tp_0_inflated_sd)*penal_tp_0 + null_penal_tp;
  nll += GMRF(QQ_tp)(tp_params);

  //likelihood for DHS data
  vector<Type> muf, mum, varf, varm;
  muf = exp(LQ_baseline_mx_DX_f + h_DX_f*h_params_f + h2_DX_f*h2_params_f + k_DX_f*k_params_f + tp_DX_f*tp_params)*Ef; 
  mum = exp(LQ_baseline_mx_DX_m + h_DX_m*h_params_m + h2_DX_m*h2_params_m + k_DX_m*k_params_m + tp_DX_m*tp_params)*Em; 

  varf = muf * (1 + muf / exp(log_dispersion(0)));
  varm = mum * (1 + mum / exp(log_dispersion(1)));
  nll -= dnbinom2(df, muf, varf, 1).sum();
  nll -= dnbinom2(dm, mum, varm, 1).sum();

  return Type(nll);

}

