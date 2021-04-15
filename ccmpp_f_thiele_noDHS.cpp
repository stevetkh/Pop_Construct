#include <TMB.hpp>                                
#include "ccmpp.h"

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
  DATA_VECTOR(log_fx_mean);
  DATA_VECTOR(srb);
  DATA_MATRIX(census_log_pop_f);
  DATA_IVECTOR(census_year_idx);
  DATA_IVECTOR(census_year_grow_idx);
  DATA_SCALAR(interval);
  DATA_INTEGER(n_periods);
  DATA_INTEGER(fx_idx);
  DATA_INTEGER(n_fx);

  DATA_INTEGER(pop_start);
  DATA_INTEGER(pop_end);
  DATA_INTEGER(open_idx);

  PARAMETER_VECTOR(log_tau2_logpop_f);
  PARAMETER(log_tau2_fx);
  PARAMETER(log_tau2_gx_f);
  PARAMETER(logit_rho_g_x);
  PARAMETER(logit_rho_g_t);

  PARAMETER_VECTOR(log_basepop_f);
  PARAMETER_VECTOR(log_fx);
  PARAMETER_VECTOR(gx_f);

  DATA_VECTOR(log_phi_mean);
  DATA_VECTOR(log_psi_mean);
  DATA_VECTOR(log_lambda_mean);
  DATA_VECTOR(log_delta_mean);
  DATA_VECTOR(log_epsilon_mean);
  DATA_VECTOR(log_A_mean);
  DATA_VECTOR(log_B_mean);

  DATA_VECTOR(thiele_age);

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

  //inverse gamma prior for variance with shape=1 and scale=0.0109
  nll -= dlgamma(log_tau2_logpop_f(0), Type(1.0), Type(1.0 / 0.0109), true);
  Type sigma_logpop_f(exp(-0.5 * log_tau2_logpop_f(0)));

  nll -= dlgamma(log_tau2_logpop_f(1), Type(1.0), Type(1.0 / 0.0109), true);
  Type sigma_logpop_f_base(exp(-0.5 * log_tau2_logpop_f(1)));

  nll -= dlgamma(log_tau2_fx, Type(1.0), Type(1.0 / 0.0109), true);
  Type sigma_fx(exp(-0.5 * log_tau2_fx));

  nll -= dlgamma(log_tau2_gx_f, Type(1.0), Type(1.0 / 0.0436), true);
  Type sigma_gx_f(exp(-0.5 * log_tau2_gx_f));

  nll -= dnorm(logit_rho_g_x, Type(0.0), Type(10.0), 1);
  nll -= dnorm(logit_rho_g_t, Type(0.0), Type(10.0), 1);
  
  nll -= dnorm(log_basepop_f, log_basepop_mean_f, sigma_logpop_f_base, true).sum();
  vector<Type> basepop_f(exp(log_basepop_f));
  
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

  Type rho_phi = 2.0 * invlogit(logit_rho_phi) - 1.0;
  Type rho_psi = 2.0 * invlogit(logit_rho_psi) - 1.0;
  Type rho_lambda = 2.0 * invlogit(logit_rho_lambda) - 1.0;
  Type rho_delta = 2.0 * invlogit(logit_rho_delta) - 1.0;
  Type rho_epsilon = 2.0 * invlogit(logit_rho_epsilon) - 1.0;
  Type rho_A = 2.0 * invlogit(logit_rho_A) - 1.0;
  Type rho_B = 2.0 * invlogit(logit_rho_B) - 1.0;
  
  Type rho_gx = 2.0 * invlogit(logit_rho_g_x) - 1.0;
  Type rho_gt = 2.0 * invlogit(logit_rho_g_t) - 1.0;

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
  
  matrix<Type> mx_mat_f(thiele_age.size(), phi.size()); 
  for(int i = 0; i < phi.size(); i++){
    mx_mat_f.col(i) = phi(i)*exp(-psi(i)*thiele_age) + lambda(i)*exp(-delta(i)*((thiele_age-epsilon(i))*(thiele_age-epsilon(i)))) + A(i)*exp(B(i)*thiele_age);
  }

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
  
  matrix<Type> mx_open_weighted_f = mx_mat_f.block(open_idx - 1, 0, mx_mat_f.rows() - open_idx + 1, n_periods).array() * weights_mat_f.array();

  matrix<Type> mx_mat_aggr_f(open_idx, n_periods);
  mx_mat_aggr_f.block(0, 0, open_idx - 1, n_periods) = mx_mat_f.block(0, 0, open_idx - 1, n_periods);
  mx_mat_aggr_f.row(open_idx - 1) = mx_open_weighted_f.colwise().sum();
  
  matrix<Type> sx_mat_f(mx_mat_aggr_f.rows() + 1, n_periods);

 for(int i=0; i < n_periods; i++) {
   //UDD 0.5q0 (2.5q[0-5] in this case)
   sx_mat_f(0,i) = Type(1.0) / (Type(1.0) + Type(0.5) * interval * mx_mat_f(0,i)); 
   //sx_mat_f(1,i) = (1.0 - exp(h_params_f(i))) / ( (1.0 - 0.5 * exp(h_params_f(i))) * (1.0 + 0.5 * interval * mx_mat_f(1,i)) ); 
 
   //UDD open age group sx
    sx_mat_f(mx_mat_aggr_f.rows(), i) = (Type(1.0) - Type(0.5) * interval * mx_mat_aggr_f(mx_mat_aggr_f.rows() - 1, i)) / ((Type(1.0) + Type(0.5) * interval * mx_mat_aggr_f(mx_mat_aggr_f.rows() - 1, i)));
    
    for(int j=0; j < (mx_mat_aggr_f.rows() - 1); j++) {
      sx_mat_f(j+1,i) = (Type(1.0) - Type(0.5) * interval * mx_mat_f(j,i)) / (Type(1.0) + Type(0.5) * interval * mx_mat_f(j+1,i));
    }
   }

  // prior for log(fx)
  nll -= dnorm(log_fx, log_fx_mean, sigma_fx, true).sum();
  vector<Type> fx(exp(log_fx));
  MapMatrixXXt fx_mat(fx.data(), n_fx, n_periods);

  // prior for gx
  MapMatrixXXt gx_mat_f(gx_f.data(), basepop_f.size(), n_periods);
  array<Type> gx_f_array(basepop_f.size(), n_periods);
  gx_f_array = gx_mat_f.array();

  nll += SCALE( SEPARABLE( AR1(rho_gt), AR1(rho_gx) ), sigma_gx_f ) (gx_f_array);

  // population projection
  PopulationProjection<Type> proj(ccmpp<Type>(basepop_f, sx_mat_f, fx_mat, gx_mat_f,
					      srb, interval, fx_idx-1));

  matrix<Type> census_proj_mat_f(basepop_f.size(),census_year_idx.size());

  for(int i = 0; i < census_year_idx.size(); i++) {	
	census_proj_mat_f.col(i) = vector<Type>(proj.population.col(census_year_idx[i] - 1)).pow(Type(1.0) - census_year_grow_idx[i] / interval) * vector<Type>(proj.population.col(census_year_idx[i])).pow(census_year_grow_idx[i] / interval);
 }

  // likelihood for log census counts
  for(int i = 0; i < census_year_idx.size(); i++) {
    nll -= dnorm(vector<Type>(census_log_pop_f.block(pop_start - 1, i, pop_end - pop_start + 1, 1)),
  		 log(vector<Type>(census_proj_mat_f.block(pop_start - 1, i, pop_end - pop_start + 1, 1))),
  		 sigma_logpop_f, true).sum();
  }

  MatrixXXt proj_period_deaths_f(proj.period_deaths());
  vector<Type> period_deaths_f(MapVectorXt(proj_period_deaths_f.data(),
					 proj_period_deaths_f.size()));
  vector<Type> births(MapVectorXt(proj.births.data(),
  				  proj.births.size()));


  DATA_INTEGER(calc_outputs);

  if(calc_outputs) {
    
    vector<Type> population_f(MapVectorXt(proj.population.data(),
					proj.population.size()));
    vector<Type> cohort_deaths_f(MapVectorXt(proj.cohort_deaths.data(),
					   proj.cohort_deaths.size()));
    vector<Type> infants_f(MapVectorXt(proj.infants.data(),
				     proj.infants.size()));
    vector<Type> migrations_f(MapVectorXt(proj.migrations.data(),
					proj.migrations.size()));
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

    REPORT(phi);
    REPORT(psi);
    REPORT(lambda);
    REPORT(delta);
    REPORT(epsilon);
    REPORT(A);
    REPORT(B);
  }

  return Type(nll);

}
