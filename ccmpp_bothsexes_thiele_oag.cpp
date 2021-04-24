#include <TMB.hpp>                                
#include "ccmpp.h"
#include "ccmpp_m.h"

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
 
  DATA_VECTOR(log_basepop_mean_f);
  DATA_VECTOR(log_basepop_mean_m);
  DATA_VECTOR(log_fx_mean);
  DATA_VECTOR(srb);
  DATA_MATRIX(census_log_pop_f);
  DATA_MATRIX(census_log_pop_m);
  DATA_IVECTOR(census_year_idx);
  DATA_IVECTOR(census_year_grow_idx);
  DATA_IVECTOR(oag);
  DATA_SCALAR(interval);
  DATA_INTEGER(n_periods);
  DATA_INTEGER(fx_idx);
  DATA_INTEGER(n_fx);
  DATA_IVECTOR(pop_start);
  DATA_IVECTOR(pop_end);
  DATA_INTEGER(open_idx);

  PARAMETER_VECTOR(log_tau2_logpop_f);
  PARAMETER_VECTOR(log_tau2_logpop_m);
  PARAMETER(log_tau2_fx);
  PARAMETER(log_tau2_gx_f);
  PARAMETER(log_tau2_gx_m);
  PARAMETER(logit_rho_g_x_f);
  PARAMETER(logit_rho_g_t_f);
  PARAMETER(logit_rho_g_x_m);
  PARAMETER(logit_rho_g_t_m);

  PARAMETER_VECTOR(log_basepop_f);
  PARAMETER_VECTOR(log_basepop_m);
  PARAMETER_VECTOR(log_fx);
  PARAMETER_VECTOR(gx_f);
  PARAMETER_VECTOR(gx_m);

  DATA_VECTOR(df);
  DATA_VECTOR(dm);
  DATA_VECTOR(Ef);
  DATA_VECTOR(Em);
  DATA_IVECTOR(df_age);
  DATA_IVECTOR(dm_age);
  DATA_IVECTOR(df_time);
  DATA_IVECTOR(dm_time);
  DATA_IVECTOR(df_tp);
  DATA_IVECTOR(dm_tp);
  DATA_VECTOR(thiele_age);

  DATA_VECTOR(log_phi_mean_f);
  DATA_VECTOR(log_psi_mean_f);
  DATA_VECTOR(log_lambda_mean_f);
  DATA_VECTOR(log_delta_mean_f);
  DATA_VECTOR(log_epsilon_mean_f);
  DATA_VECTOR(log_A_mean_f);
  DATA_VECTOR(log_B_mean_f);

  DATA_VECTOR(log_phi_mean_m);
  DATA_VECTOR(log_psi_mean_m);
  DATA_VECTOR(log_lambda_mean_m);
  DATA_VECTOR(log_delta_mean_m);
  DATA_VECTOR(log_epsilon_mean_m);
  DATA_VECTOR(log_A_mean_m);
  DATA_VECTOR(log_B_mean_m);

  DATA_SPARSE_MATRIX(penal_tp);
  DATA_SPARSE_MATRIX(penal_tp_0);
  DATA_SPARSE_MATRIX(null_penal_tp);

  PARAMETER(log_lambda_tp);
  PARAMETER(log_lambda_tp_0_inflated_sd);
  PARAMETER(log_dispersion_f);
  PARAMETER(log_dispersion_m);

  PARAMETER_VECTOR(tp_params);

  PARAMETER_VECTOR(log_phi_innov_f); 
  PARAMETER_VECTOR(log_psi_innov_f);
  PARAMETER_VECTOR(log_lambda_innov_f);
  PARAMETER_VECTOR(log_delta_innov_f);
  PARAMETER_VECTOR(log_epsilon_innov_f);
  PARAMETER_VECTOR(log_A_innov_f);
  PARAMETER_VECTOR(log_B_innov_f);

  PARAMETER_VECTOR(log_phi_innov_m); 
  PARAMETER_VECTOR(log_psi_innov_m);
  PARAMETER_VECTOR(log_lambda_innov_m);
  PARAMETER_VECTOR(log_delta_innov_m);
  PARAMETER_VECTOR(log_epsilon_innov_m);
  PARAMETER_VECTOR(log_A_innov_m);
  PARAMETER_VECTOR(log_B_innov_m);

  PARAMETER(log_marginal_prec_phi_f);
  PARAMETER(log_marginal_prec_psi_f);
  PARAMETER(log_marginal_prec_lambda_f);
  PARAMETER(log_marginal_prec_delta_f);
  PARAMETER(log_marginal_prec_epsilon_f);
  PARAMETER(log_marginal_prec_A_f);
  PARAMETER(log_marginal_prec_B_f);

  PARAMETER(log_marginal_prec_phi_m);
  PARAMETER(log_marginal_prec_psi_m);
  PARAMETER(log_marginal_prec_lambda_m);
  PARAMETER(log_marginal_prec_delta_m);
  PARAMETER(log_marginal_prec_epsilon_m);
  PARAMETER(log_marginal_prec_A_m);
  PARAMETER(log_marginal_prec_B_m);

  PARAMETER(logit_rho_phi_f);
  PARAMETER(logit_rho_psi_f);
  PARAMETER(logit_rho_lambda_f);
  PARAMETER(logit_rho_delta_f);
  PARAMETER(logit_rho_epsilon_f);
  PARAMETER(logit_rho_A_f);
  PARAMETER(logit_rho_B_f);

  PARAMETER(logit_rho_phi_m);
  PARAMETER(logit_rho_psi_m);
  PARAMETER(logit_rho_lambda_m);
  PARAMETER(logit_rho_delta_m);
  PARAMETER(logit_rho_epsilon_m);
  PARAMETER(logit_rho_A_m);
  PARAMETER(logit_rho_B_m);

  Type nll(0.0);

  //inverse gamma prior for variance with shape=1 and scale=0.0109  
  nll -= dlgamma(log_tau2_logpop_f(0), Type(1.0), Type(1.0 / 0.0109), true);
  Type sigma_logpop_f(exp(-0.5 * log_tau2_logpop_f(0)));
  nll -= dlgamma(log_tau2_logpop_m(0), Type(1.0), Type(1.0 / 0.0109), true);
  Type sigma_logpop_m(exp(-0.5 * log_tau2_logpop_m(0)));

  nll -= dlgamma(log_tau2_logpop_f(1), Type(1.0), Type(1.0 / 0.0109), true);
  Type sigma_logpop_f_base(exp(-0.5 * log_tau2_logpop_f(1)));
  nll -= dlgamma(log_tau2_logpop_m(1), Type(1.0), Type(1.0 / 0.0109), true);
  Type sigma_logpop_m_base(exp(-0.5 * log_tau2_logpop_m(1)));

  nll -= dlgamma(log_tau2_fx, Type(1.0), Type(1.0 / 0.0109), true);
  Type sigma_fx(exp(-0.5 * log_tau2_fx));

  nll -= dlgamma(log_tau2_gx_f, Type(1.0), Type(1.0 / 0.0436), true);
  Type sigma_gx_f(exp(-0.5 * log_tau2_gx_f));
  nll -= dlgamma(log_tau2_gx_m, Type(1.0), Type(1.0 / 0.0436), true);
  Type sigma_gx_m(exp(-0.5 * log_tau2_gx_m));

  nll -= dnorm(logit_rho_g_x_f, Type(0.0), Type(5.0), 1);
  nll -= dnorm(logit_rho_g_t_f, Type(0.0), Type(5.0), 1);
  nll -= dnorm(logit_rho_g_x_m, Type(0.0), Type(5.0), 1);
  nll -= dnorm(logit_rho_g_t_m, Type(0.0), Type(5.0), 1);  

  nll -= dnorm(log_basepop_f, log_basepop_mean_f, sigma_logpop_f_base, true).sum();
  vector<Type> basepop_f(exp(log_basepop_f));
  nll -= dnorm(log_basepop_m, log_basepop_mean_m, sigma_logpop_m_base, true).sum();
  vector<Type> basepop_m(exp(log_basepop_m));

  nll -= dlgamma(log_marginal_prec_phi_f, Type(1.0), Type(1.0 / 0.01), true);
  Type sigma_phi_f(exp(-0.5 * log_marginal_prec_phi_f));
  nll -= dlgamma(log_marginal_prec_phi_m, Type(1.0), Type(1.0 / 0.01), true);
  Type sigma_phi_m(exp(-0.5 * log_marginal_prec_phi_m));
  
  nll -= dlgamma(log_marginal_prec_psi_f, Type(1.0), Type(1.0 / 0.01), true);
  Type sigma_psi_f(exp(-0.5 * log_marginal_prec_psi_f));
  nll -= dlgamma(log_marginal_prec_psi_m, Type(1.0), Type(1.0 / 0.01), true);
  Type sigma_psi_m(exp(-0.5 * log_marginal_prec_psi_m));  

  nll -= dlgamma(log_marginal_prec_lambda_f, Type(1.0), Type(1.0 / 0.01), true);
  Type sigma_lambda_f(exp(-0.5 * log_marginal_prec_lambda_f));
  nll -= dlgamma(log_marginal_prec_lambda_m, Type(1.0), Type(1.0 / 0.01), true);
  Type sigma_lambda_m(exp(-0.5 * log_marginal_prec_lambda_m));
  
  nll -= dlgamma(log_marginal_prec_delta_f, Type(1.0), Type(1.0 / 0.01), true);
  Type sigma_delta_f(exp(-0.5 * log_marginal_prec_delta_f));
  nll -= dlgamma(log_marginal_prec_delta_m, Type(1.0), Type(1.0 / 0.01), true);
  Type sigma_delta_m(exp(-0.5 * log_marginal_prec_delta_m));  

  nll -= dlgamma(log_marginal_prec_epsilon_f, Type(1.0), Type(1.0 / 0.01), true);
  Type sigma_epsilon_f(exp(-0.5 * log_marginal_prec_epsilon_f));
  nll -= dlgamma(log_marginal_prec_epsilon_m, Type(1.0), Type(1.0 / 0.01), true);
  Type sigma_epsilon_m(exp(-0.5 * log_marginal_prec_epsilon_m));
  
  nll -= dlgamma(log_marginal_prec_A_f, Type(1.0), Type(1.0 / 0.01), true);
  Type sigma_A_f(exp(-0.5 * log_marginal_prec_A_f));
  nll -= dlgamma(log_marginal_prec_A_m, Type(1.0), Type(1.0 / 0.01), true);
  Type sigma_A_m(exp(-0.5 * log_marginal_prec_A_m));

  nll -= dlgamma(log_marginal_prec_B_f, Type(1.0), Type(1.0 / 0.01), true);
  Type sigma_B_f(exp(-0.5 * log_marginal_prec_B_f));
  nll -= dlgamma(log_marginal_prec_B_m, Type(1.0), Type(1.0 / 0.01), true);
  Type sigma_B_m(exp(-0.5 * log_marginal_prec_B_m));

  nll -= dnorm(logit_rho_phi_f, Type(0.0), Type(10.0), 1);
  nll -= dnorm(logit_rho_psi_f, Type(0.0), Type(10.0), 1);
  nll -= dnorm(logit_rho_lambda_f, Type(0.0), Type(10.0), 1);
  nll -= dnorm(logit_rho_delta_f, Type(0.0), Type(10.0), 1);
  nll -= dnorm(logit_rho_epsilon_f, Type(0.0), Type(10.0), 1);
  nll -= dnorm(logit_rho_A_f, Type(0.0), Type(10.0), 1);
  nll -= dnorm(logit_rho_B_f, Type(0.0), Type(10.0), 1);

  nll -= dnorm(logit_rho_phi_m, Type(0.0), Type(10.0), 1);
  nll -= dnorm(logit_rho_psi_m, Type(0.0), Type(10.0), 1);
  nll -= dnorm(logit_rho_lambda_m, Type(0.0), Type(10.0), 1);
  nll -= dnorm(logit_rho_delta_m, Type(0.0), Type(10.0), 1);
  nll -= dnorm(logit_rho_epsilon_m, Type(0.0), Type(10.0), 1);
  nll -= dnorm(logit_rho_A_m, Type(0.0), Type(10.0), 1);
  nll -= dnorm(logit_rho_B_m, Type(0.0), Type(10.0), 1);

  nll -= dnorm(log_lambda_tp, Type(0.0), Type(5.0), 1);
  nll -= dnorm(log_lambda_tp_0_inflated_sd, Type(0.0), Type(5.0), 1);

  nll -= dnorm(log_dispersion_f, Type(0.0), Type(5.0), 1);
  nll -= dnorm(log_dispersion_m, Type(0.0), Type(5.0), 1);

  nll -= dnorm(log_phi_innov_f, Type(0.0), Type(1.0), 1).sum();
  nll -= dnorm(log_psi_innov_f, Type(0.0), Type(1.0), 1).sum();
  nll -= dnorm(log_lambda_innov_f, Type(0.0), Type(1.0), 1).sum();
  nll -= dnorm(log_delta_innov_f, Type(0.0), Type(1.0), 1).sum();
  nll -= dnorm(log_epsilon_innov_f, Type(0.0), Type(1.0), 1).sum();
  nll -= dnorm(log_A_innov_f, Type(0.0), Type(1.0), 1).sum();
  nll -= dnorm(log_B_innov_f, Type(0.0), Type(1.0), 1).sum();

  nll -= dnorm(log_phi_innov_m, Type(0.0), Type(1.0), 1).sum();
  nll -= dnorm(log_psi_innov_m, Type(0.0), Type(1.0), 1).sum();
  nll -= dnorm(log_lambda_innov_m, Type(0.0), Type(1.0), 1).sum();
  nll -= dnorm(log_delta_innov_m, Type(0.0), Type(1.0), 1).sum();
  nll -= dnorm(log_epsilon_innov_m, Type(0.0), Type(1.0), 1).sum();
  nll -= dnorm(log_A_innov_m, Type(0.0), Type(1.0), 1).sum();
  nll -= dnorm(log_B_innov_m, Type(0.0), Type(1.0), 1).sum();

  Type rho_phi_f = 2.0 * invlogit(logit_rho_phi_f) - 1.0;
  Type rho_psi_f = 2.0 * invlogit(logit_rho_psi_f) - 1.0;
  Type rho_lambda_f = 2.0 * invlogit(logit_rho_lambda_f) - 1.0;
  Type rho_delta_f = 2.0 * invlogit(logit_rho_delta_f) - 1.0;
  Type rho_epsilon_f = 2.0 * invlogit(logit_rho_epsilon_f) - 1.0;
  Type rho_A_f = 2.0 * invlogit(logit_rho_A_f) - 1.0;
  Type rho_B_f = 2.0 * invlogit(logit_rho_B_f) - 1.0;

  Type rho_phi_m = 2.0 * invlogit(logit_rho_phi_m) - 1.0;
  Type rho_psi_m = 2.0 * invlogit(logit_rho_psi_m) - 1.0;
  Type rho_lambda_m = 2.0 * invlogit(logit_rho_lambda_m) - 1.0;
  Type rho_delta_m = 2.0 * invlogit(logit_rho_delta_m) - 1.0;
  Type rho_epsilon_m = 2.0 * invlogit(logit_rho_epsilon_m) - 1.0;
  Type rho_A_m = 2.0 * invlogit(logit_rho_A_m) - 1.0;
  Type rho_B_m = 2.0 * invlogit(logit_rho_B_m) - 1.0;

  Type rho_gx_f = 2.0 * invlogit(logit_rho_g_x_f) - 1.0;
  Type rho_gt_f = 2.0 * invlogit(logit_rho_g_t_f) - 1.0;
  Type rho_gx_m = 2.0 * invlogit(logit_rho_g_x_m) - 1.0;
  Type rho_gt_m = 2.0 * invlogit(logit_rho_g_t_m) - 1.0;

  vector<Type> log_phi_f(log_phi_innov_f.size());
  vector<Type> log_psi_f(log_phi_innov_f.size());
  vector<Type> log_lambda_f(log_phi_innov_f.size());
  vector<Type> log_delta_f(log_phi_innov_f.size());
  vector<Type> log_epsilon_f(log_phi_innov_f.size());
  vector<Type> log_A_f(log_phi_innov_f.size());
  vector<Type> log_B_f(log_phi_innov_f.size());

  vector<Type> log_phi_m(log_phi_innov_m.size());
  vector<Type> log_psi_m(log_phi_innov_m.size());
  vector<Type> log_lambda_m(log_phi_innov_m.size());
  vector<Type> log_delta_m(log_phi_innov_m.size());
  vector<Type> log_epsilon_m(log_phi_innov_m.size());
  vector<Type> log_A_m(log_phi_innov_m.size());
  vector<Type> log_B_m(log_phi_innov_m.size());

  log_phi_f(0) = sigma_phi_f * log_phi_innov_f(0);
  log_psi_f(0) = sigma_psi_f * log_psi_innov_f(0);
  log_lambda_f(0) = sigma_lambda_f * log_lambda_innov_f(0);
  log_delta_f(0) = sigma_delta_f * log_delta_innov_f(0);
  log_epsilon_f(0) = sigma_epsilon_f * log_epsilon_innov_f(0);
  log_A_f(0) = sigma_A_f * log_A_innov_f(0);
  log_B_f(0) = sigma_B_f * log_B_innov_f(0);

  log_phi_m(0) = sigma_phi_m * log_phi_innov_m(0);
  log_psi_m(0) = sigma_psi_m * log_psi_innov_m(0);
  log_lambda_m(0) = sigma_lambda_m * log_lambda_innov_m(0);
  log_delta_m(0) = sigma_delta_m * log_delta_innov_m(0);
  log_epsilon_m(0) = sigma_epsilon_m * log_epsilon_innov_m(0);
  log_A_m(0) = sigma_A_m * log_A_innov_m(0);
  log_B_m(0) = sigma_B_m * log_B_innov_m(0);

  for(int i=1; i < log_phi_innov_f.size(); i++){
    log_phi_f(i) = rho_phi_f * log_phi_f(i-1) + sqrt(1.0 - rho_phi_f * rho_phi_f) * sigma_phi_f * log_phi_innov_f(i);
    log_psi_f(i) = rho_psi_f * log_psi_f(i-1) + sqrt(1.0 - rho_psi_f * rho_psi_f) * sigma_psi_f * log_psi_innov_f(i);
    log_lambda_f(i) = rho_lambda_f * log_lambda_f(i-1) + sqrt(1.0 - rho_lambda_f * rho_lambda_f) * sigma_lambda_f * log_lambda_innov_f(i);
    log_delta_f(i) = rho_delta_f * log_delta_f(i-1) + sqrt(1.0 - rho_delta_f * rho_delta_f) * sigma_delta_f * log_delta_innov_f(i);
    log_epsilon_f(i) = rho_epsilon_f * log_epsilon_f(i-1) + sqrt(1.0 - rho_epsilon_f * rho_epsilon_f) * sigma_epsilon_f * log_epsilon_innov_f(i);
    log_A_f(i) = rho_A_f * log_A_f(i-1) + sqrt(1.0 - rho_A_f * rho_A_f) * sigma_A_f * log_A_innov_f(i);
    log_B_f(i) = rho_B_f * log_B_f(i-1) + sqrt(1.0 - rho_B_f * rho_B_f) * sigma_B_f * log_B_innov_f(i);

    log_phi_m(i) = rho_phi_m * log_phi_m(i-1) + sqrt(1.0 - rho_phi_m * rho_phi_m) * sigma_phi_m * log_phi_innov_m(i);
    log_psi_m(i) = rho_psi_m * log_psi_m(i-1) + sqrt(1.0 - rho_psi_m * rho_psi_m) * sigma_psi_m * log_psi_innov_m(i);
    log_lambda_m(i) = rho_lambda_m * log_lambda_m(i-1) + sqrt(1.0 - rho_lambda_m * rho_lambda_m) * sigma_lambda_m * log_lambda_innov_m(i);
    log_delta_m(i) = rho_delta_m * log_delta_m(i-1) + sqrt(1.0 - rho_delta_m * rho_delta_m) * sigma_delta_m * log_delta_innov_m(i);
    log_epsilon_m(i) = rho_epsilon_m * log_epsilon_m(i-1) + sqrt(1.0 - rho_epsilon_m * rho_epsilon_m) * sigma_epsilon_m * log_epsilon_innov_m(i);
    log_A_m(i) = rho_A_m * log_A_m(i-1) + sqrt(1.0 - rho_A_m * rho_A_m) * sigma_A_m * log_A_innov_m(i);
    log_B_m(i) = rho_B_m * log_B_m(i-1) + sqrt(1.0 - rho_B_m * rho_B_m) * sigma_B_m * log_B_innov_m(i);
  }

  log_phi_f += log_phi_mean_f;  
  log_psi_f += log_psi_mean_f;  
  log_lambda_f += log_lambda_mean_f;  
  log_delta_f += log_delta_mean_f;  
  log_epsilon_f += log_epsilon_mean_f;  
  log_A_f += log_A_mean_f;  
  log_B_f += log_B_mean_f;  

  log_phi_m += log_phi_mean_m;  
  log_psi_m += log_psi_mean_m;  
  log_lambda_m += log_lambda_mean_m;  
  log_delta_m += log_delta_mean_m;  
  log_epsilon_m += log_epsilon_mean_m;  
  log_A_m += log_A_mean_m;  
  log_B_m += log_B_mean_m;  

  vector<Type> phi_f = exp(log_phi_f);
  vector<Type> psi_f = exp(log_psi_f);
  vector<Type> lambda_f = exp(log_lambda_f);
  vector<Type> delta_f = exp(log_delta_f);
  vector<Type> epsilon_f = exp(log_epsilon_f);
  vector<Type> A_f = exp(log_A_f);
  vector<Type> B_f = exp(log_B_f);

  vector<Type> phi_m = exp(log_phi_m);
  vector<Type> psi_m = exp(log_psi_m);
  vector<Type> lambda_m = exp(log_lambda_m);
  vector<Type> delta_m = exp(log_delta_m);
  vector<Type> epsilon_m = exp(log_epsilon_m);
  vector<Type> A_m = exp(log_A_m);
  vector<Type> B_m = exp(log_B_m);

  matrix<Type> mx_mat_f(thiele_age.size(), phi_f.size()); 
  matrix<Type> mx_mat_m(thiele_age.size(), phi_f.size()); 
  for(int i = 0; i < phi_f.size(); i++){
    mx_mat_f.col(i) = phi_f(i)*exp(-psi_f(i)*(thiele_age - Type(2.0))) + lambda_f(i)*exp(-delta_f(i)*((thiele_age-epsilon_f(i))*(thiele_age-epsilon_f(i)))) + A_f(i)*exp(B_f(i)*(thiele_age - Type(92.0)));
    mx_mat_m.col(i) = phi_m(i)*exp(-psi_m(i)*(thiele_age - Type(2.0))) + lambda_m(i)*exp(-delta_m(i)*((thiele_age-epsilon_m(i))*(thiele_age-epsilon_m(i)))) + A_m(i)*exp(B_m(i)*(thiele_age - Type(92.0)));
  }

  SparseMatrix<Type> QQ_tp = exp(log_lambda_tp)*penal_tp + exp(-2*log_lambda_tp_0_inflated_sd)*penal_tp_0 + null_penal_tp;
  nll += GMRF(QQ_tp)(tp_params);

  //likelihood for DHS data
  vector<Type> muf(df.size());
  vector<Type> varf(df.size());
  vector<Type> mum(dm.size());
  vector<Type> varm(dm.size());
  for(int i = 0; i < df.size(); i++){
    muf(i) = mx_mat_f(df_age(i)-1, df_time(i)-1) * exp(tp_params(df_tp(i))) * Ef(i);
  }
  for(int i = 0; i < dm.size(); i++){
    mum(i) = mx_mat_m(dm_age(i)-1, dm_time(i)-1) * exp(tp_params(dm_tp(i))) * Em(i);
  }

  varf = muf * (1 + muf / exp(log_dispersion_f));
  varm = mum * (1 + mum / exp(log_dispersion_m));
  nll -= dnbinom2(df, muf, varf, 1).sum();
  nll -= dnbinom2(dm, mum, varm, 1).sum();

  matrix<Type> weights_mat_f(mx_mat_f.rows() - open_idx + 1, n_periods);
  weights_mat_f = Type(1.0) / (Type(1.0) + Type(0.5) * interval * mx_mat_f.block(open_idx - 1, 0, mx_mat_f.rows() - open_idx + 1, n_periods).array());
  
  for(int i=0; i < n_periods; i++) {
  	for(int j=open_idx - 1; j < mx_mat_f.rows()-1; j++) {
		for(int k=j-open_idx+2; k < mx_mat_f.rows()-open_idx+1; k++) {
  			weights_mat_f(k,i) *= (Type(1.0) - Type(0.5) * interval * mx_mat_f(j,i)) / (Type(1.0) + Type(0.5) * interval * mx_mat_f(j,i));
		}
  	}
  }
  vector<Type> weights_sum_f = weights_mat_f.colwise().sum();
  for(int i=0; i < n_periods; i++) {
    for(int j=0; j < weights_mat_f.rows(); j++) {
    weights_mat_f(j,i) *= 1 / weights_sum_f(i);
    }
  }
  
  matrix<Type> weights_mat_m(mx_mat_m.rows() - open_idx + 1, n_periods);
  weights_mat_m = Type(1.0) / (Type(1.0) + Type(0.5) * interval * mx_mat_m.block(open_idx - 1, 0, mx_mat_m.rows() - open_idx + 1, n_periods).array());
  
  for(int i=0; i < n_periods; i++) {
  	for(int j=open_idx - 1; j < mx_mat_m.rows()-1; j++) {
		for(int k=j-open_idx+2; k < mx_mat_m.rows()-open_idx+1; k++) {
  			weights_mat_m(k,i) *= (Type(1.0) - Type(0.5) * interval * mx_mat_m(j,i)) / (Type(1.0) + Type(0.5) * interval * mx_mat_m(j,i));
		}
  	}
  }
  vector<Type> weights_sum_m = weights_mat_m.colwise().sum();
  for(int i=0; i < n_periods; i++) {
    for(int j=0; j < weights_mat_m.rows(); j++) {
    weights_mat_m(j,i) *= 1 / weights_sum_m(i);
    }
  }
  

  matrix<Type> mx_open_weighted_f = mx_mat_f.block(open_idx - 1, 0, mx_mat_f.rows() - open_idx + 1, n_periods).array() * weights_mat_f.array();
  matrix<Type> mx_open_weighted_m = mx_mat_m.block(open_idx - 1, 0, mx_mat_m.rows() - open_idx + 1, n_periods).array() * weights_mat_m.array();

  matrix<Type> mx_mat_aggr_f(open_idx, n_periods);
  mx_mat_aggr_f.block(0, 0, open_idx - 1, n_periods) = mx_mat_f.block(0, 0, open_idx - 1, n_periods);
  mx_mat_aggr_f.row(open_idx - 1) = mx_open_weighted_f.colwise().sum();
  
  matrix<Type> mx_mat_aggr_m(open_idx, n_periods);
  mx_mat_aggr_m.block(0, 0, open_idx - 1, n_periods) = mx_mat_m.block(0, 0, open_idx - 1, n_periods);
  mx_mat_aggr_m.row(open_idx - 1) = mx_open_weighted_m.colwise().sum();
  
  matrix<Type> sx_mat_f(mx_mat_aggr_f.rows() + 1, n_periods);
  matrix<Type> sx_mat_m(mx_mat_aggr_m.rows() + 1, n_periods);
  
  for(int i=0; i < n_periods; i++) {
    //UDD 0.5q0 (2.5q[0-5] in this case)
    sx_mat_f(0,i) = Type(1.0) / (Type(1.0) + Type(0.5) * interval * mx_mat_f(0,i));
    sx_mat_m(0,i) = Type(1.0) / (Type(1.0) + Type(0.5) * interval * mx_mat_m(0,i));
    //UDD open age group sx
    sx_mat_f(mx_mat_aggr_f.rows(), i) = (1 - 2.5*mx_mat_aggr_f(mx_mat_aggr_f.rows() - 1, i)) / ((1 + 2.5*mx_mat_aggr_f(mx_mat_aggr_f.rows() - 1, i)));
    sx_mat_m(mx_mat_aggr_m.rows(), i) = (1 - 2.5*mx_mat_aggr_m(mx_mat_aggr_m.rows() - 1, i)) / ((1 + 2.5*mx_mat_aggr_m(mx_mat_aggr_m.rows() - 1, i)));
    
    for(int j=0; j < (mx_mat_aggr_f.rows() - 1); j++) {
      sx_mat_f(j+1,i) = (Type(1.0) - Type(0.5) * interval * mx_mat_f(j,i)) / (Type(1.0) + Type(0.5) * interval * mx_mat_f(j+1,i));
      sx_mat_m(j+1,i) = (Type(1.0) - Type(0.5) * interval * mx_mat_m(j,i)) / (Type(1.0) + Type(0.5) * interval * mx_mat_m(j+1,i));
     }
   }

  // prior for log(fx)
  nll -= dnorm(log_fx, log_fx_mean, sigma_fx, true).sum();
  vector<Type> fx(exp(log_fx));
  MapMatrixXXt fx_mat(fx.data(), n_fx, n_periods);

  // prior for gx
  MapMatrixXXt gx_mat_f(gx_f.data(), basepop_f.size(), n_periods);
  MapMatrixXXt gx_mat_m(gx_m.data(), basepop_m.size(), n_periods);
  array<Type> gx_f_array(basepop_f.size(), n_periods);
  array<Type> gx_m_array(basepop_m.size(), n_periods);
  gx_f_array = gx_mat_f.array();
  gx_m_array = gx_mat_m.array();

  nll += SCALE( SEPARABLE( AR1(rho_gt_f), AR1(rho_gx_f) ), sigma_gx_f ) (gx_f_array);
  nll += SCALE( SEPARABLE( AR1(rho_gt_m), AR1(rho_gx_m) ), sigma_gx_m ) (gx_m_array);

  // population projection
  PopulationProjection<Type> proj(ccmpp<Type>(basepop_f, sx_mat_f, fx_mat, gx_mat_f,
					      srb, interval, fx_idx-1));
  PopulationProjection_m<Type> proj_m(ccmpp_m<Type>(basepop_m, sx_mat_m, gx_mat_m,
					      srb, interval, proj.births));

  matrix<Type> census_proj_mat_f(basepop_f.size(),census_year_idx.size());
  matrix<Type> census_proj_mat_m(basepop_m.size(),census_year_idx.size());
 
  matrix<Type> census_proj_mat_f_oag(oag.maxCoeff(),census_year_idx.size());
  matrix<Type> census_proj_mat_m_oag(oag.maxCoeff(),census_year_idx.size());

  for(int i = 0; i < census_year_idx.size(); i++) {	
	census_proj_mat_f.col(i) = vector<Type>(proj.population.col(census_year_idx[i] - 1)).pow(Type(1.0) - census_year_grow_idx[i] / interval) * vector<Type>(proj.population.col(census_year_idx[i])).pow(census_year_grow_idx[i] / interval);
	census_proj_mat_m.col(i) = vector<Type>(proj_m.population.col(census_year_idx[i] - 1)).pow(Type(1.0) - census_year_grow_idx[i] / interval) * vector<Type>(proj_m.population.col(census_year_idx[i])).pow(census_year_grow_idx[i] / interval);

	census_proj_mat_f_oag.block(0, i, oag(i)-1, 1) = census_proj_mat_f.block(0, i, oag(i)-1, 1);
	census_proj_mat_m_oag.block(0, i, oag(i)-1, 1) = census_proj_mat_m.block(0, i, oag(i)-1, 1);
	census_proj_mat_f_oag(oag(i)-1, i)  = census_proj_mat_f.col(i).segment(oag(i)-1, basepop_f.size() - oag(i) + 1).sum();
	census_proj_mat_m_oag(oag(i)-1, i)  = census_proj_mat_m.col(i).segment(oag(i)-1, basepop_m.size() - oag(i) + 1).sum();
  }
 
  // likelihood for log census counts
  for(int i = 0; i < census_year_idx.size(); i++) {
    nll -= dnorm(vector<Type>(census_log_pop_f.block(pop_start(i) - 1, i, pop_end(i) - pop_start(i) + 1, 1)),
  		 log(vector<Type>(census_proj_mat_f_oag.block(pop_start(i) - 1, i, pop_end(i) - pop_start(i) + 1, 1))),
  		 sigma_logpop_f, true).sum();
    nll -= dnorm(vector<Type>(census_log_pop_m.block(pop_start(i) - 1, i, pop_end(i) - pop_start(i) + 1, 1)),
  		 log(vector<Type>(census_proj_mat_m_oag.block(pop_start(i) - 1, i, pop_end(i) - pop_start(i) + 1, 1))),
  		 sigma_logpop_m, true).sum();
  }

  MatrixXXt proj_period_deaths_f(proj.period_deaths());
  MatrixXXt proj_period_deaths_m(proj_m.period_deaths());
  vector<Type> period_deaths_f(MapVectorXt(proj_period_deaths_f.data(),
					 proj_period_deaths_f.size()));
  vector<Type> period_deaths_m(MapVectorXt(proj_period_deaths_m.data(),
					 proj_period_deaths_m.size()));
  vector<Type> births(MapVectorXt(proj.births.data(),
  				  proj.births.size()));


  DATA_INTEGER(calc_outputs);

  if(calc_outputs) {
    
    vector<Type> population_f(MapVectorXt(proj.population.data(),
					proj.population.size()));
    vector<Type> population_m(MapVectorXt(proj_m.population.data(),
					proj_m.population.size()));
    vector<Type> cohort_deaths_f(MapVectorXt(proj.cohort_deaths.data(),
					   proj.cohort_deaths.size()));
    vector<Type> cohort_deaths_m(MapVectorXt(proj_m.cohort_deaths.data(),
					   proj_m.cohort_deaths.size()));
    vector<Type> infants_f(MapVectorXt(proj.infants.data(),
				     proj.infants.size()));
    vector<Type> infants_m(MapVectorXt(proj_m.infants.data(),
				     proj_m.infants.size()));
    vector<Type> migrations_f(MapVectorXt(proj.migrations.data(),
					proj.migrations.size()));
    vector<Type> migrations_m(MapVectorXt(proj_m.migrations.data(),
					proj_m.migrations.size()));
    REPORT(population_f);
    REPORT(cohort_deaths_f);
    REPORT(period_deaths_f);
    REPORT(births);
    REPORT(infants_f);
    REPORT(migrations_f);
    REPORT(sx_mat_f);
    REPORT(mx_mat_f);
    REPORT(fx);
    REPORT(gx_f);
    REPORT(census_proj_mat_f);

    REPORT(population_m);
    REPORT(cohort_deaths_m);
    REPORT(period_deaths_m);
    REPORT(infants_m);
    REPORT(migrations_m);
    REPORT(sx_mat_m);
    REPORT(mx_mat_m);
    REPORT(gx_m);
    REPORT(census_proj_mat_m);

    REPORT(phi_f);
    REPORT(psi_f);
    REPORT(lambda_f);
    REPORT(delta_f);
    REPORT(epsilon_f);
    REPORT(A_f);
    REPORT(B_f);

    REPORT(phi_m);
    REPORT(psi_m);
    REPORT(lambda_m);
    REPORT(delta_m);
    REPORT(epsilon_m);
    REPORT(A_m);
    REPORT(B_m);
  }

  return Type(nll);

}
