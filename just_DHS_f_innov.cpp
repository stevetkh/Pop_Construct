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

  PARAMETER(log_marginal_prec_h);
  PARAMETER(log_marginal_prec_k);
  PARAMETER(log_lambda_tp);
  PARAMETER(log_lambda_tp_0_inflated_sd);
  PARAMETER(log_dispersion);
  PARAMETER(logit_rho_h);
  PARAMETER(logit_rho_k);
 
  PARAMETER_VECTOR(h_params_f_innov);
  PARAMETER_VECTOR(k_params_f_innov);

  PARAMETER_VECTOR(h_constant_f);
  PARAMETER_VECTOR(tp_params);

  Type nll(0.0);

  nll -= dlgamma(log_marginal_prec_h, Type(1.0), Type(1.0 / 0.0109), true);
  Type sigma_h(exp(-0.5 * log_marginal_prec_h));
  nll -= dlgamma(log_marginal_prec_k, Type(1.0), Type(1.0 / 2), true);
  Type sigma_k(exp(-0.5 * log_marginal_prec_k));

//  nll -= dnorm(log_marginal_prec_h, Type(0.0), Type(5.0), true);
//  Type sigma_h(exp(0.5 * log_marginal_prec_h));

//  nll -= dnorm(log_marginal_prec_k, Type(0.0), Type(5.0), true);
//  Type sigma_k(exp(0.5 * log_marginal_prec_k));

  nll -= dnorm(logit_rho_h, Type(0.0), Type(10.0), 1);
  nll -= dnorm(logit_rho_k, Type(0.0), Type(10.0), 1);

  nll -= dnorm(log_lambda_tp, Type(0.0), Type(5.0), 1);
  nll -= dnorm(log_lambda_tp_0_inflated_sd, Type(0.0), Type(5.0), 1);
  
  nll -= dnorm(log_dispersion, Type(0.0), Type(5.0), 1);

  Type rho_h = 2.0 * invlogit(logit_rho_h) - 1.0;
  Type rho_k = 2.0 * invlogit(logit_rho_k) - 1.0;
   
//  nll -= dnorm(rho_h, Type(0.0), Type(0.5), 1);
//  nll -= dnorm(rho_k, Type(0.0), Type(0.5), 1);
//  nll -= logit_rho_h - 2 * log(1 + exp(logit_rho_h));
//  nll -= logit_rho_k - 2 * log(1 + exp(logit_rho_k));

  nll -= dnorm(h_params_f_innov, Type(0.0), Type(1.0), 1).sum();
  nll -= dnorm(k_params_f_innov, Type(0.0), Type(1.0), 1).sum();

  vector<Type> h_params_f(h_params_f_innov.size());
  vector<Type> k_params_f(k_params_f_innov.size());

  h_params_f(0) = sigma_h * h_params_f_innov(0);
  k_params_f(0) = sigma_k * k_params_f_innov(0);

  for(int i=1; i < h_params_f_innov.size(); i++){
    h_params_f(i) = rho_h * h_params_f(i-1) + sqrt(1.0 - rho_h * rho_h) * sigma_h * h_params_f_innov(i);
    k_params_f(i) = rho_k * k_params_f(i-1) + sqrt(1.0 - rho_k * rho_k) * sigma_k * k_params_f_innov(i);
  }

  h_params_f += h_constant_f;

  SparseMatrix<Type> QQ_tp = exp(log_lambda_tp)*penal_tp + exp(-2*log_lambda_tp_0_inflated_sd)*penal_tp_0 + null_penal_tp;
  nll += GMRF(QQ_tp)(tp_params);

  vector<Type> h2_params_f = h_params_f*h_params_f;

  //likelihood for DHS data
  vector<Type> muf, varf;
  muf = exp(LQ_baseline_mx_DX_f + h_DX_f*h_params_f + h2_DX_f*h2_params_f + k_DX_f*k_params_f + tp_DX_f*tp_params)*Ef; 

  varf = muf * (1 + muf / exp(log_dispersion));
  nll -= dnbinom2(df, muf, varf, 1).sum();

  REPORT(h_params_f);
  REPORT(k_params_f);
  	
  return Type(nll);

}
