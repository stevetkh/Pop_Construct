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
  DATA_VECTOR(DEF_log_age);
  DATA_VECTOR(h_constant_f);
  DATA_SPARSE_MATRIX(DEF_time);

  PARAMETER(log_marginal_prec_h);
//  PARAMETER(log_marginal_prec_k);
  PARAMETER(log_lambda_tp);
  PARAMETER(log_lambda_tp_0_inflated_sd);
  PARAMETER(log_dispersion);
  PARAMETER(logit_rho_h);
 
  PARAMETER_VECTOR(h_params_f_innov);
  PARAMETER(k_param_f);

  PARAMETER_VECTOR(tp_params);

//  PARAMETER(log_D);
//  PARAMETER(log_E);
//  PARAMETER(log_F);

  PARAMETER_VECTOR(log_D_innov);
  PARAMETER(log_D_mean);
  PARAMETER_VECTOR(log_E_innov);
  PARAMETER(log_E_mean);
  PARAMETER_VECTOR(log_F_innov);
  PARAMETER(log_F_mean);

  PARAMETER(logit_rho_D);
  PARAMETER(logit_rho_E);
  PARAMETER(logit_rho_F);

  PARAMETER(log_marginal_prec_D);
  PARAMETER(log_marginal_prec_E);
  PARAMETER(log_marginal_prec_F);

  Type nll(0.0);

  nll -= dlgamma(log_marginal_prec_h, Type(1.0), Type(1.0 / 0.0109), true);
  Type sigma_h(exp(-0.5 * log_marginal_prec_h));
  
//  nll -= dlgamma(log_marginal_prec_k, Type(1.0), Type(1.0 / 2.0), true);
//  Type sigma_k(exp(-0.5 * log_marginal_prec_k));

  nll -= dlgamma(log_marginal_prec_D, Type(1.0), Type(1.0 / 5.0), true);
  Type sigma_D(exp(-0.5 * log_marginal_prec_D));

  nll -= dlgamma(log_marginal_prec_E, Type(1.0), Type(1.0 / 5.0), true);
  Type sigma_E(exp(-0.5 * log_marginal_prec_E));

  nll -= dlgamma(log_marginal_prec_F, Type(1.0), Type(1.0 / 5.0), true);
  Type sigma_F(exp(-0.5 * log_marginal_prec_F));

  nll -= dnorm(logit_rho_h, Type(0.0), Type(10.0), 1);

  nll -= dnorm(logit_rho_D, Type(0.0), Type(10.0), 1);
  nll -= dnorm(logit_rho_E, Type(0.0), Type(10.0), 1);
  nll -= dnorm(logit_rho_F, Type(0.0), Type(10.0), 1);

//  nll -= dnorm(log_D_mean, Type(-6.0), sigma_D_mean, 1);
//  nll -= dnorm(log_E_mean, Type(1.5), sigma_E_mean, 1);
//  nll -= dnorm(log_F_mean, Type(3.5), sigma_F_mean, 1);

  nll -= dnorm(log_lambda_tp, Type(0.0), Type(5.0), 1);
  nll -= dnorm(log_lambda_tp_0_inflated_sd, Type(0.0), Type(5.0), 1);
  
  nll -= dnorm(log_dispersion, Type(0.0), Type(5.0), 1);

  Type rho_h = 2.0 * invlogit(logit_rho_h) - 1.0;
  Type rho_D = 2.0 * invlogit(logit_rho_D) - 1.0;
  Type rho_E = 2.0 * invlogit(logit_rho_E) - 1.0;
  Type rho_F = 2.0 * invlogit(logit_rho_F) - 1.0;
  
  nll -= dnorm(h_params_f_innov, Type(0.0), Type(1.0), 1).sum();
  nll -= dnorm(k_param_f, Type(0.0), Type(1.0), 1);

  nll -= dnorm(log_D_innov, Type(0.0), Type(1.0), 1).sum();
  nll -= dnorm(log_E_innov, Type(0.0), Type(1.0), 1).sum();
  nll -= dnorm(log_F_innov, Type(0.0), Type(1.0), 1).sum();

  vector<Type> h_params_f(h_params_f_innov.size());
  vector<Type> k_params_f(h_params_f_innov.size());
  vector<Type> log_D(h_params_f_innov.size());
  vector<Type> log_E(h_params_f_innov.size());
  vector<Type> log_F(h_params_f_innov.size());

  h_params_f(0) = sigma_h * h_params_f_innov(0);
  k_params_f(0) = k_param_f;
  log_D(0) = sigma_D * log_D_innov(0);
  log_E(0) = sigma_E * log_E_innov(0);
  log_F(0) = sigma_F * log_F_innov(0);

  for(int i=1; i < h_params_f_innov.size(); i++){
    h_params_f(i) = rho_h * h_params_f(i-1) + sqrt(1.0 - rho_h * rho_h) * sigma_h * h_params_f_innov(i);
    k_params_f(i) = k_param_f;
    log_D(i) = rho_D * log_D(i-1) + sqrt(1.0 - rho_D * rho_D) * sigma_D * log_D_innov(i);
    log_E(i) = rho_E * log_E(i-1) + sqrt(1.0 - rho_E * rho_E) * sigma_E * log_E_innov(i);
    log_F(i) = rho_F * log_F(i-1) + sqrt(1.0 - rho_F * rho_F) * sigma_F * log_F_innov(i);
  }

  h_params_f += h_constant_f;
  vector<Type> h2_params_f = h_params_f*h_params_f;
  
  log_D += log_D_mean;
  log_E += log_E_mean;
  log_F += log_F_mean;

  SparseMatrix<Type> QQ_tp = exp(log_lambda_tp)*penal_tp + exp(-2*log_lambda_tp_0_inflated_sd)*penal_tp_0 + null_penal_tp;
  nll += GMRF(QQ_tp)(tp_params);

  vector<Type> log_hump;
  log_hump = DEF_log_age - DEF_time * log_F;
  log_hump *= DEF_log_age - DEF_time * log_F;	
  log_hump *= -exp(DEF_time * log_E);
  
  //likelihood for DHS data
  vector<Type> muf, varf;
  muf = exp(LQ_baseline_mx_DX_f + h_DX_f*h_params_f + h2_DX_f*h2_params_f + k_DX_f*k_params_f + tp_DX_f*tp_params)*Ef; 
  muf += exp(DEF_time * log_D + log_hump + tp_DX_f*tp_params)*Ef;

  varf = muf * (1 + muf / exp(log_dispersion));
  nll -= dnbinom2(df, muf, varf, 1).sum();

  REPORT(h_params_f);
  REPORT(k_params_f);
  REPORT(log_D);
  REPORT(log_E);
  REPORT(log_F);

  return Type(nll);

}
