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
 
  DATA_VECTOR(df);
  DATA_VECTOR(Ef);
  DATA_IVECTOR(df_age);
  DATA_IVECTOR(df_time);
  DATA_IVECTOR(df_tp);

  DATA_VECTOR(log_phi_mean);
  DATA_VECTOR(log_psi_mean);
  DATA_VECTOR(log_lambda_mean);
  DATA_VECTOR(log_delta_mean);
  DATA_VECTOR(log_epsilon_mean);
  DATA_VECTOR(log_A_mean);
  DATA_VECTOR(log_B_mean);

  DATA_VECTOR(thiele_age);

  DATA_SPARSE_MATRIX(penal_tp);
  DATA_SPARSE_MATRIX(penal_tp_0);
  DATA_SPARSE_MATRIX(null_penal_tp);

  PARAMETER(log_lambda_tp);
  PARAMETER(log_lambda_tp_0_inflated_sd);
  PARAMETER_VECTOR(tp_params);

  PARAMETER(log_dispersion);

  PARAMETER_VECTOR(log_phi_innov); 
  PARAMETER_VECTOR(log_psi_innov);
  PARAMETER_VECTOR(log_lambda_innov);
  PARAMETER_VECTOR(log_delta_innov);
  PARAMETER_VECTOR(log_epsilon_innov);
  PARAMETER_VECTOR(log_A_innov);
  PARAMETER_VECTOR(log_B_innov);

  PARAMETER(log_marginal_prec_phi);
  PARAMETER(log_marginal_prec_psi);
  PARAMETER(log_marginal_prec_lambda);
  PARAMETER(log_marginal_prec_delta);
  PARAMETER(log_marginal_prec_epsilon);
  PARAMETER(log_marginal_prec_A);
  PARAMETER(log_marginal_prec_B);

  PARAMETER(logit_rho_phi);
  PARAMETER(logit_rho_psi);
  PARAMETER(logit_rho_lambda);
  PARAMETER(logit_rho_delta);
  PARAMETER(logit_rho_epsilon);
  PARAMETER(logit_rho_A);
  PARAMETER(logit_rho_B);

  Type nll(0.0);

  nll -= dlgamma(log_marginal_prec_phi, Type(1.0), Type(1.0 / 0.01), true);
  Type sigma_phi(exp(-0.5 * log_marginal_prec_phi));
  
  nll -= dlgamma(log_marginal_prec_psi, Type(1.0), Type(1.0 / 0.01), true);
  Type sigma_psi(exp(-0.5 * log_marginal_prec_psi));
  
  nll -= dlgamma(log_marginal_prec_lambda, Type(1.0), Type(1.0 / 0.01), true);
  Type sigma_lambda(exp(-0.5 * log_marginal_prec_lambda));
  
  nll -= dlgamma(log_marginal_prec_delta, Type(1.0), Type(1.0 / 0.01), true);
  Type sigma_delta(exp(-0.5 * log_marginal_prec_delta));
  
  nll -= dlgamma(log_marginal_prec_epsilon, Type(1.0), Type(1.0 / 0.01), true);
  Type sigma_epsilon(exp(-0.5 * log_marginal_prec_epsilon));
  
  nll -= dlgamma(log_marginal_prec_A, Type(1.0), Type(1.0 / 0.01), true);
  Type sigma_A(exp(-0.5 * log_marginal_prec_A));

  nll -= dlgamma(log_marginal_prec_B, Type(1.0), Type(1.0 / 0.01), true);
  Type sigma_B(exp(-0.5 * log_marginal_prec_B));
  
  nll -= dnorm(logit_rho_phi, Type(0.0), Type(10.0), 1);
  nll -= dnorm(logit_rho_psi, Type(0.0), Type(10.0), 1);
  nll -= dnorm(logit_rho_lambda, Type(0.0), Type(10.0), 1);
  nll -= dnorm(logit_rho_delta, Type(0.0), Type(10.0), 1);
  nll -= dnorm(logit_rho_epsilon, Type(0.0), Type(10.0), 1);
  nll -= dnorm(logit_rho_A, Type(0.0), Type(10.0), 1);
  nll -= dnorm(logit_rho_B, Type(0.0), Type(10.0), 1);

  nll -= dnorm(log_lambda_tp, Type(0.0), Type(5.0), 1);
  nll -= dnorm(log_lambda_tp_0_inflated_sd, Type(0.0), Type(5.0), 1);
  
  nll -= dnorm(log_dispersion, Type(0.0), Type(5.0), 1);

  Type rho_phi = 2.0 * invlogit(logit_rho_phi) - 1.0;
  Type rho_psi = 2.0 * invlogit(logit_rho_psi) - 1.0;
  Type rho_lambda = 2.0 * invlogit(logit_rho_lambda) - 1.0;
  Type rho_delta = 2.0 * invlogit(logit_rho_delta) - 1.0;
  Type rho_epsilon = 2.0 * invlogit(logit_rho_epsilon) - 1.0;
  Type rho_A = 2.0 * invlogit(logit_rho_A) - 1.0;
  Type rho_B = 2.0 * invlogit(logit_rho_B) - 1.0;
  
  nll -= dnorm(log_phi_innov, Type(0.0), Type(1.0), 1).sum();
  nll -= dnorm(log_psi_innov, Type(0.0), Type(1.0), 1).sum();
  nll -= dnorm(log_lambda_innov, Type(0.0), Type(1.0), 1).sum();
  nll -= dnorm(log_delta_innov, Type(0.0), Type(1.0), 1).sum();
  nll -= dnorm(log_epsilon_innov, Type(0.0), Type(1.0), 1).sum();
  nll -= dnorm(log_A_innov, Type(0.0), Type(1.0), 1).sum();
  nll -= dnorm(log_B_innov, Type(0.0), Type(1.0), 1).sum();

  vector<Type> log_phi(log_phi_innov.size());
  vector<Type> log_psi(log_phi_innov.size());
  vector<Type> log_lambda(log_phi_innov.size());
  vector<Type> log_delta(log_phi_innov.size());
  vector<Type> log_epsilon(log_phi_innov.size());
  vector<Type> log_A(log_phi_innov.size());
  vector<Type> log_B(log_phi_innov.size());

  log_phi(0) = sigma_phi * log_phi_innov(0);
  log_psi(0) = sigma_psi * log_psi_innov(0);
  log_lambda(0) = sigma_lambda * log_lambda_innov(0);
  log_delta(0) = sigma_delta * log_delta_innov(0);
  log_epsilon(0) = sigma_epsilon * log_epsilon_innov(0);
  log_A(0) = sigma_A * log_A_innov(0);
  log_B(0) = sigma_B * log_B_innov(0);

  for(int i=1; i < log_phi_innov.size(); i++){
    log_phi(i) = rho_phi * log_phi(i-1) + sqrt(1.0 - rho_phi * rho_phi) * sigma_phi * log_phi_innov(i);
    log_psi(i) = rho_psi * log_psi(i-1) + sqrt(1.0 - rho_psi * rho_psi) * sigma_psi * log_psi_innov(i);
    log_lambda(i) = rho_lambda * log_lambda(i-1) + sqrt(1.0 - rho_lambda * rho_lambda) * sigma_lambda * log_lambda_innov(i);
    log_delta(i) = rho_delta * log_delta(i-1) + sqrt(1.0 - rho_delta * rho_delta) * sigma_delta * log_delta_innov(i);
    log_epsilon(i) = rho_epsilon * log_epsilon(i-1) + sqrt(1.0 - rho_epsilon * rho_epsilon) * sigma_epsilon * log_epsilon_innov(i);
    log_A(i) = rho_A * log_A(i-1) + sqrt(1.0 - rho_A * rho_A) * sigma_A * log_A_innov(i);
    log_B(i) = rho_B * log_B(i-1) + sqrt(1.0 - rho_B * rho_B) * sigma_B * log_B_innov(i);
  }

  log_phi += log_phi_mean;  
  log_psi += log_psi_mean;  
  log_lambda += log_lambda_mean;  
  log_delta += log_delta_mean;  
  log_epsilon += log_epsilon_mean;  
  log_A += log_A_mean;  
  log_B += log_B_mean;  

  vector<Type> phi = exp(log_phi);
  vector<Type> psi = exp(log_psi);
  vector<Type> lambda = exp(log_lambda);
  vector<Type> delta = exp(log_delta);
  vector<Type> epsilon = exp(log_epsilon);
  vector<Type> A = exp(log_A);
  vector<Type> B = exp(log_B);
  
  matrix<Type> mx_mat(thiele_age.size(), phi.size()); 
  for(int i = 0; i < phi.size(); i++){
    mx_mat.col(i) = phi(i)*exp(-psi(i)*thiele_age) + lambda(i)*exp(-delta(i)*((thiele_age-epsilon(i))*(thiele_age-epsilon(i)))) + A(i)*exp(B(i)*thiele_age);
  }

  SparseMatrix<Type> QQ_tp = exp(log_lambda_tp)*penal_tp + exp(-2*log_lambda_tp_0_inflated_sd)*penal_tp_0 + null_penal_tp;
  nll += GMRF(QQ_tp)(tp_params);

  //likelihood for DHS data
  vector<Type> muf(df.size());
  vector<Type> varf(df.size());
  for(int i = 0; i < df.size(); i++){
    muf(i) = mx_mat(df_age(i)-1, df_time(i)-1) * exp(tp_params(df_tp(i))) * Ef(i);
  }
 
  varf = muf * (1 + muf / exp(log_dispersion));
  nll -= dnbinom2(df, muf, varf, 1).sum();

  REPORT(phi);
  REPORT(psi);
  REPORT(lambda);
  REPORT(delta);
  REPORT(epsilon);
  REPORT(A);
  REPORT(B);
  REPORT(mx_mat);

  return Type(nll);

}

