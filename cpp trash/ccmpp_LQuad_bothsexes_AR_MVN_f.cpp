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

  DATA_INTEGER(open_idx);
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
  DATA_INTEGER(pop_start);
  DATA_INTEGER(pop_end);

  PARAMETER(log_tau2_logpop_f);
  PARAMETER(log_tau2_fx);
  PARAMETER(log_tau2_gx_f);
  PARAMETER(log_marginal_prec_h);
  PARAMETER(log_marginal_prec_k);
  PARAMETER(log_lambda_tp);
  PARAMETER(log_lambda_tp_0_inflated_sd);
  PARAMETER(log_dispersion);
 
  PARAMETER(logit_rho_g_x);
  PARAMETER(logit_rho_g_t);

  PARAMETER_VECTOR(log_basepop_f);
  PARAMETER_VECTOR(log_fx);
  PARAMETER_VECTOR(gx_f);
  PARAMETER_VECTOR(h_params_f);
  PARAMETER_VECTOR(k_params_f);
  PARAMETER(logit_rho_h);
  PARAMETER_VECTOR(h_constant_f);
  PARAMETER_VECTOR(tp_params);
  vector<Type> h2_params_f = h_params_f*h_params_f;

  Type nll(0.0);

  //inverse gamma prior for variance with shape=1 and scale=0.0109
  nll -= dlgamma(log_tau2_logpop_f, Type(1.0), Type(1.0 / 0.0109), true);
  Type sigma_logpop_f(exp(-0.5 * log_tau2_logpop_f));

  nll -= dlgamma(log_tau2_fx, Type(1.0), Type(1.0 / 0.0109), true);
  Type sigma_fx(exp(-0.5 * log_tau2_fx));

  nll -= dlgamma(log_tau2_gx_f, Type(1.0), Type(1.0 / 0.0436), true);
  Type sigma_gx_f(exp(-0.5 * log_tau2_gx_f));

  nll -= dnorm(logit_rho_g_x, Type(0.0), Type(10.0), 1);
  nll -= dnorm(logit_rho_g_t, Type(0.0), Type(10.0), 1);
  
  nll -= dlgamma(log_marginal_prec_h, Type(1.0), Type(1.0 / 0.0109), true);
  Type sigma_h(exp(-0.5 * log_marginal_prec_h));
  
  nll -= dlgamma(log_marginal_prec_k, Type(1.0), Type(1.0 / 2), true);
  Type sigma_k(exp(-0.5 * log_marginal_prec_k));

  nll -= dnorm(logit_rho_h, Type(0.0), Type(10.0), 1);

  nll -= dnorm(h_constant_f, Type(-2.0), Type(10.0), 1).sum();

  nll -= dnorm(log_lambda_tp, Type(0.0), Type(10.0), 1);
  nll -= dnorm(log_lambda_tp_0_inflated_sd, Type(0.0), Type(10.0), 1);
  
  nll -= dnorm(log_dispersion, Type(0.0), Type(10.0), 1);

  nll -= dnorm(log_basepop_f, log_basepop_mean_f, sigma_logpop_f, true).sum();
  vector<Type> basepop_f(exp(log_basepop_f));

  nll += SCALE(AR1(invlogit(logit_rho_h)), sigma_h)(h_params_f - h_constant_f);
  nll -= dnorm(k_params_f, Type(0.0), sigma_k, 1).sum();

  SparseMatrix<Type> QQ_tp = exp(log_lambda_tp)*penal_tp + exp(-2*log_lambda_tp_0_inflated_sd)*penal_tp_0 + null_penal_tp;
  nll += GMRF(QQ_tp)(tp_params);

  vector<Type> log_mort_rates_0 = h_params_f - log(1 - 0.5 * exp(h_params_f));

  //likelihood for DHS data
  vector<Type> muf, varf;
  muf = exp(LQ_baseline_mx_DX_f + h_DX_f*h_params_f + h2_DX_f*h2_params_f + k_DX_f*k_params_f + tp_DX_f*tp_params)*Ef; 

  varf = muf * (1 + muf / exp(log_dispersion));
  nll -= dnbinom2(df, muf, varf, 1).sum();

  matrix<Type> LQ_params_mat_f(4, n_periods);
  for(int i=0; i < n_periods; i++){
  LQ_params_mat_f(0,i) = Type(1.0);
  }
  LQ_params_mat_f.row(1) = h_params_f;
  LQ_params_mat_f.row(2) = h2_params_f;
  LQ_params_mat_f.row(3) = k_params_f;

  matrix<Type> log_mx_mat_f_no0 = LQ_baseline_f*LQ_params_mat_f;
  matrix<Type> log_mx_mat_f(log_mx_mat_f_no0.rows() + 1, log_mx_mat_f_no0.cols());
  log_mx_mat_f.bottomRows(log_mx_mat_f_no0.rows()) = log_mx_mat_f_no0;
  log_mx_mat_f.row(0) = log_mort_rates_0;

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
  
  matrix<Type> mx_open_weighted_f = mx_mat_f.block(open_idx - 1, 0, mx_mat_f.rows() - open_idx + 1, n_periods).array() * weights_mat_f.array();

  matrix<Type> mx_mat_aggr_f(open_idx, n_periods);
  mx_mat_aggr_f.block(0, 0, open_idx - 1, n_periods) = mx_mat_f.block(0, 0, open_idx - 1, n_periods);
  mx_mat_aggr_f.row(open_idx - 1) = mx_open_weighted_f.colwise().sum();
  
  matrix<Type> sx_mat_f(mx_mat_aggr_f.rows() + 1, n_periods);

  for(int i=0; i < n_periods; i++) {
    //UDD 0.5q0 (2.5q[0-5] in this case)
    sx_mat_f(0,i) = Type(1.0) / (Type(1.0) + Type(0.5) * interval * mx_mat_aggr_f(0,i)); 
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

  nll += SCALE( SEPARABLE( AR1(invlogit(logit_rho_g_t)), AR1(invlogit(logit_rho_g_x)) ), sigma_gx_f ) (gx_f_array);

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
  }

  return Type(nll);

}
