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
  DATA_VECTOR(log_fx_mean);
  DATA_VECTOR(srb);
  DATA_MATRIX(census_log_pop_f);
  DATA_IVECTOR(census_year_idx);
  DATA_IVECTOR(census_year_grow_idx);
  DATA_SCALAR(interval);
  DATA_INTEGER(n_periods);
  DATA_INTEGER(fx_idx);
  DATA_INTEGER(n_fx);

  DATA_VECTOR(log_basepop_mean_m);
  DATA_MATRIX(census_log_pop_m);

  DATA_INTEGER(open_idx);
  DATA_INTEGER(pop_start);
  DATA_INTEGER(pop_end); 
  DATA_MATRIX(LQ_baseline_f);
  DATA_MATRIX(LQ_baseline_m);

  PARAMETER(log_tau2_logpop_f);
  PARAMETER(log_tau2_fx);
  PARAMETER(log_tau2_gx_f);
  PARAMETER_VECTOR(log_marginal_prec_h);
  PARAMETER_VECTOR(log_marginal_prec_k);
 
  PARAMETER(log_tau2_logpop_m);
  PARAMETER(log_tau2_gx_m);

  PARAMETER_VECTOR(logit_rho_g_x);
  PARAMETER_VECTOR(logit_rho_g_t);

  PARAMETER_VECTOR(log_basepop_f);
  PARAMETER_VECTOR(log_fx);
  PARAMETER_VECTOR(gx_f);
  PARAMETER_VECTOR(h_params_f);
  PARAMETER_VECTOR(k_params_f);
  PARAMETER_VECTOR(logit_rho_h);
  PARAMETER_VECTOR(logit_rho_k);
  PARAMETER_VECTOR(h_constant_f);
  PARAMETER(k_constant_f);
  vector<Type> h2_params_f = h_params_f*h_params_f;

  PARAMETER_VECTOR(log_basepop_m);
  PARAMETER_VECTOR(gx_m);
  PARAMETER_VECTOR(h_params_m);
  PARAMETER_VECTOR(k_params_m);
  PARAMETER_VECTOR(h_constant_m);
  PARAMETER(k_constant_m);
  vector<Type> h2_params_m = h_params_m*h_params_m;

  Type nll(0.0);

  nll -= dlgamma(log_tau2_logpop_f, Type(1.0), Type(1.0 / 0.0109), true);
  Type sigma_logpop_f(exp(-0.5 * log_tau2_logpop_f));
  nll -= dlgamma(log_tau2_logpop_m, Type(1.0), Type(1.0 / 0.0109), true);
  Type sigma_logpop_m(exp(-0.5 * log_tau2_logpop_m));

  nll -= dlgamma(log_tau2_fx, Type(1.0), Type(1.0 / 0.0109), true);
  Type sigma_fx(exp(-0.5 * log_tau2_fx));

  nll -= dlgamma(log_tau2_gx_f, Type(1.0), Type(1.0 / 0.0436), true);
  Type sigma_gx_f(exp(-0.5 * log_tau2_gx_f));
  nll -= dlgamma(log_tau2_gx_m, Type(1.0), Type(1.0 / 0.0436), true);
  Type sigma_gx_m(exp(-0.5 * log_tau2_gx_m));

  nll -= dnorm(logit_rho_g_x, Type(0.0), Type(10.0), 1).sum();
  nll -= dnorm(logit_rho_g_t, Type(0.0), Type(10.0), 1).sum();

  nll -= dlgamma(log_marginal_prec_h, Type(1.0), Type(1.0 / 0.0109), true).sum();
  vector<Type> sigma_h(exp(-0.5 * log_marginal_prec_h));

  nll -= dlgamma(log_marginal_prec_k, Type(1.0), Type(1.0 / 2), true).sum();
  vector<Type> sigma_k(exp(-0.5 * log_marginal_prec_k));

  nll -= dnorm(logit_rho_h, Type(0.0), Type(10.0), 1).sum();
  nll -= dnorm(logit_rho_k, Type(0.0), Type(10.0), 1).sum();

  nll -= dnorm(h_constant_f, Type(-2.0), Type(10.0), 1).sum();
  nll -= dnorm(h_constant_m, Type(-2.0), Type(10.0), 1).sum();

  nll -= dnorm(k_constant_f, Type(0.0), Type(10.0), 1);
  nll -= dnorm(k_constant_m, Type(0.0), Type(10.0), 1);

  nll -= dnorm(log_basepop_f, log_basepop_mean_f, sigma_logpop_f, true).sum();
  vector<Type> basepop_f(exp(log_basepop_f));
  nll -= dnorm(log_basepop_m, log_basepop_mean_m, sigma_logpop_m, true).sum();
  vector<Type> basepop_m(exp(log_basepop_m));

  vector<Type> log_mort_rates_0_f = h_params_f - log(1.0 - 0.5 * exp(h_params_f)) + log(1.0 / interval);
  vector<Type> log_mort_rates_0_m = h_params_m - log(1.0 - 0.5 * exp(h_params_m)) + log(1.0 / interval);

  vector<Type> rho_h = 2.0 * invlogit(logit_rho_h) - 1.0;
  vector<Type> rho_k = 2.0 * invlogit(logit_rho_k) - 1.0;
  vector<Type> rho_gx = 2.0 * invlogit(logit_rho_g_x) - 1.0;
  vector<Type> rho_gt = 2.0 * invlogit(logit_rho_g_t) - 1.0;

  nll += SCALE(AR1(rho_h(0)), sigma_h(0))(h_params_f - h_constant_f);
  nll += SCALE(AR1(rho_h(1)), sigma_h(1))(h_params_m - h_constant_m);
  nll += SCALE(AR1(rho_k(0)), sigma_k(0))(k_params_f - k_constant_f);
  nll += SCALE(AR1(rho_k(1)), sigma_k(1))(k_params_m - k_constant_m);

  matrix<Type> LQ_params_mat_f(4, n_periods);
  matrix<Type> LQ_params_mat_m(4, n_periods);
  for(int i=0; i < n_periods; i++){
  LQ_params_mat_f(0,i) = Type(1.0);
  LQ_params_mat_m(0,i) = Type(1.0);
  }
  LQ_params_mat_f.row(1) = h_params_f;
  LQ_params_mat_f.row(2) = h2_params_f;
  LQ_params_mat_f.row(3) = k_params_f;
  LQ_params_mat_m.row(1) = h_params_m;
  LQ_params_mat_m.row(2) = h2_params_m;
  LQ_params_mat_m.row(3) = k_params_m;

  matrix<Type> log_mx_mat_f(LQ_baseline_f.rows() + 1, n_periods);
  log_mx_mat_f.row(0) = log_mort_rates_0_f;
  log_mx_mat_f.bottomRows(LQ_baseline_f.rows()) = LQ_baseline_f*LQ_params_mat_f;

  matrix<Type> mx_mat_f = exp(log_mx_mat_f.array());
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
  
  matrix<Type> log_mx_mat_m(LQ_baseline_m.rows() + 1, n_periods);
  log_mx_mat_m.row(0) = log_mort_rates_0_m;
  log_mx_mat_m.bottomRows(LQ_baseline_m.rows()) = LQ_baseline_m*LQ_params_mat_m;

  matrix<Type> mx_mat_m = exp(log_mx_mat_m.array());
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
    weights_mat_m(j,i) *= Type(1.0) / weights_sum_m(i);
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
   sx_mat_f(0,i) = 1.0 - 0.5 * exp(h_params_f(i)); 
   sx_mat_f(1,i) = (1.0 - exp(h_params_f(i))) / ( (1.0 - 0.5 * exp(h_params_f(i))) * (1.0 + 0.5 * interval * mx_mat_f(1,i)) ); 
   sx_mat_m(0,i) = 1.0 - 0.5 * exp(h_params_m(i)); 
   sx_mat_m(1,i) = (1.0 - exp(h_params_m(i))) / ( (1.0 - 0.5 * exp(h_params_m(i))) * (1.0 + 0.5 * interval * mx_mat_m(1,i)) ); 
 
    //UDD open age group sx
    sx_mat_f(mx_mat_aggr_f.rows(), i) = (Type(1.0) - Type(0.5) * interval * mx_mat_aggr_f(mx_mat_aggr_f.rows() - 1, i)) / ((Type(1.0) + Type(0.5) * interval * mx_mat_aggr_f(mx_mat_aggr_f.rows() - 1, i)));
    sx_mat_m(mx_mat_aggr_m.rows(), i) = (Type(1.0) - Type(0.5) * interval * mx_mat_aggr_m(mx_mat_aggr_m.rows() - 1, i)) / ((Type(1.0) + Type(0.5) * interval * mx_mat_aggr_m(mx_mat_aggr_m.rows() - 1, i)));
    
    for(int j=1; j < (mx_mat_aggr_f.rows() - 1); j++) {
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
  gx_f_array = gx_mat_f.array();
  array<Type> gx_m_array(basepop_m.size(), n_periods);
  gx_m_array = gx_mat_m.array();

  nll += SCALE( SEPARABLE( AR1(rho_gt(0)), AR1(rho_gx(0)) ), sigma_gx_f ) (gx_f_array);
  nll += SCALE( SEPARABLE( AR1(rho_gt(1)), AR1(rho_gx(1)) ), sigma_gx_m ) (gx_m_array);

  // population projection
  PopulationProjection<Type> proj(ccmpp<Type>(basepop_f, sx_mat_f, fx_mat, gx_mat_f,
					      srb, interval, fx_idx-1));
  PopulationProjection_m<Type> proj_m(ccmpp_m<Type>(basepop_m, sx_mat_m, gx_mat_m,
					      srb, interval, proj.births));

  matrix<Type> census_proj_mat_f(basepop_f.size(),census_year_idx.size());
  matrix<Type> census_proj_mat_m(basepop_m.size(),census_year_idx.size());

  for(int i = 0; i < census_year_idx.size(); i++) {	
	census_proj_mat_f.col(i) = vector<Type>(proj.population.col(census_year_idx[i] - 1)).pow(Type(1.0) - census_year_grow_idx[i] / interval) * vector<Type>(proj.population.col(census_year_idx[i])).pow(census_year_grow_idx[i] / interval);
  	census_proj_mat_m.col(i) = vector<Type>(proj_m.population.col(census_year_idx[i] - 1)).pow(Type(1.0) - census_year_grow_idx[i] / interval) * vector<Type>(proj_m.population.col(census_year_idx[i])).pow(census_year_grow_idx[i] / interval);
 }



  // likelihood for log census counts
  for(int i = 0; i < census_year_idx.size(); i++) {
    nll -= dnorm(vector<Type>(census_log_pop_f.block(pop_start - 1, i, pop_end - pop_start + 1, 1)),
  		 log(vector<Type>(census_proj_mat_f.block(pop_start - 1, i, pop_end - pop_start + 1, 1))),
  		 sigma_logpop_f, true).sum();
    nll -= dnorm(vector<Type>(census_log_pop_m.block(pop_start - 1, i, pop_end - pop_start + 1, 1)),
  		 log(vector<Type>(census_proj_mat_m.block(pop_start - 1, i, pop_end - pop_start + 1, 1))),
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
  }

  return Type(nll);

}
