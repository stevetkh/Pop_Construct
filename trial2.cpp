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

  PARAMETER_VECTOR(log_phi_innov); 
  PARAMETER_VECTOR(log_psi_innov);
  PARAMETER_VECTOR(log_lambda_innov);
  PARAMETER_VECTOR(log_delta_innov);
  PARAMETER_VECTOR(log_epsilon_innov);
  PARAMETER_VECTOR(log_A_innov);
  PARAMETER_VECTOR(log_B_innov);

  PARAMETER(logit_rho_phi);  

  Type nll(0.0);

  nll -= dnorm(log_phi_innov, Type(0.0), Type(1.0), 1).sum();
  nll -= dnorm(log_psi_innov, Type(0.0), Type(1.0), 1).sum();
  nll -= dnorm(log_lambda_innov, Type(0.0), Type(1.0), 1).sum();
  nll -= dnorm(log_delta_innov, Type(0.0), Type(1.0), 1).sum();
  nll -= dnorm(log_epsilon_innov, Type(0.0), Type(1.0), 1).sum();
  nll -= dnorm(log_A_innov, Type(0.0), Type(1.0), 1).sum();
  nll -= dnorm(log_B_innov, Type(0.0), Type(1.0), 1).sum();

  nll -= dnorm(logit_rho_phi, Type(0.0), Type(5.0), 1);

  vector<Type> log_phi(log_phi_innov.size());
  vector<Type> log_psi(log_phi_innov.size());
  vector<Type> log_lambda(log_phi_innov.size());
  vector<Type> log_delta(log_phi_innov.size());
  vector<Type> log_epsilon(log_phi_innov.size());
  vector<Type> log_A(log_phi_innov.size());
  vector<Type> log_B(log_phi_innov.size());
  
  log_phi = log_phi_innov;
  log_psi = log_psi_innov;
  log_lambda = log_lambda_innov;
  log_delta = log_delta_innov;
  log_epsilon = log_epsilon_innov;
  log_A = log_A_innov;
  log_B = log_B_innov;
  
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

  vector<Type> muf(df.size());
  for(int i = 0; i < df.size(); i++){
    muf(i) = mx_mat(df_age(i)-1, df_time(i)-1) * Ef(i);
  }

  return Type(nll);
}


