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
 
  DATA_VECTOR(log_basepop_mean);
  DATA_VECTOR(log_fx_mean);
  DATA_VECTOR(gx_mean);
  DATA_VECTOR(srb);
  DATA_MATRIX(census_log_pop);
  DATA_IVECTOR(census_year_idx);
  DATA_SCALAR(interval);
  DATA_INTEGER(n_periods);
  DATA_INTEGER(fx_idx);
  DATA_INTEGER(n_fx);

  DATA_INTEGER(open_idx);
  DATA_VECTOR(LQ_baseline_mx_DX);
  DATA_SPARSE_MATRIX(h_DX);
  DATA_SPARSE_MATRIX(h2_DX);
  DATA_SPARSE_MATRIX(k_DX);
  DATA_SPARSE_MATRIX(tp_DX);
  DATA_SPARSE_MATRIX(penal_tp);
  DATA_SPARSE_MATRIX(penal_tp_0);
  DATA_SPARSE_MATRIX(null_penal_tp);
  DATA_SPARSE_MATRIX(penal_mx);
  DATA_SPARSE_MATRIX(penal_mx_0);
  DATA_MATRIX(LQ_baseline);
  DATA_VECTOR(dx);
  DATA_VECTOR(Ex);

  PARAMETER(log_tau2_logpop);
  PARAMETER(log_tau2_fx);
  PARAMETER(log_tau2_gx);
  PARAMETER(log_lambda_h);
  PARAMETER(log_lambda_h_0_inflated_sd);
  PARAMETER(log_lambda_k);
  PARAMETER(log_lambda_k_0_inflated_sd);
  PARAMETER(log_lambda_tp);
  PARAMETER(log_lambda_tp_0_inflated_sd);
  PARAMETER(log_dispersion);
 
  PARAMETER_VECTOR(log_basepop);
  PARAMETER_VECTOR(log_fx);
  PARAMETER_VECTOR(gx);
  PARAMETER_VECTOR(h_params);
  PARAMETER_VECTOR(k_params);
  PARAMETER_VECTOR(tp_params);
  vector<Type> h2_params = h_params*h_params;

  Type nll(0.0);

  nll -= dlgamma(log_tau2_logpop, Type(1.0), Type(1.0 / 0.0109), true);
  Type sigma_logpop(exp(-0.5 * log_tau2_logpop));

  nll -= dlgamma(log_tau2_fx, Type(1.0), Type(1.0 / 0.0109), true);
  Type sigma_fx(exp(-0.5 * log_tau2_fx));

  nll -= dlgamma(log_tau2_gx, Type(1.0), Type(1.0 / 0.0436), true);
  Type sigma_gx(exp(-0.5 * log_tau2_gx));

  nll -= dnorm(log_lambda_h, Type(0.0), Type(5.0), 1);
  nll -= dnorm(log_lambda_h_0_inflated_sd, Type(0.0), Type(5.0), 1);
 
  nll -= dnorm(log_lambda_k, Type(0.0), Type(5.0), 1);
  nll -= dnorm(log_lambda_k_0_inflated_sd, Type(0.0), Type(5.0), 1);

  nll -= dnorm(log_lambda_tp, Type(0.0), Type(5.0), 1);
  nll -= dnorm(log_lambda_tp_0_inflated_sd, Type(0.0), Type(5.0), 1);
  
  nll -= dnorm(log_dispersion, Type(0.0), Type(5.0), 1);

  nll -= dnorm(log_basepop, log_basepop_mean, sigma_logpop, true).sum();
  vector<Type> basepop(exp(log_basepop));

  //AR(1) on mx_params and tips
  SparseMatrix<Type> QQ_h = exp(log_lambda_h)*penal_mx + exp(-2*log_lambda_h_0_inflated_sd)*penal_mx_0;
  nll += GMRF(QQ_h)(h_params);
  SparseMatrix<Type> QQ_k = exp(log_lambda_k)*penal_mx + exp(-2*log_lambda_k_0_inflated_sd)*penal_mx_0;
  nll += GMRF(QQ_k)(k_params);
  SparseMatrix<Type> QQ_tp = exp(log_lambda_tp)*penal_tp + exp(-2*log_lambda_tp_0_inflated_sd)*penal_tp_0 + null_penal_tp;
  nll += GMRF(QQ_tp)(tp_params);

  //likelihood for DHS data
  vector<Type> mu, var;
  mu = exp(LQ_baseline_mx_DX + h_DX*h_params + h2_DX*h2_params + k_DX*k_params + tp_DX*tp_params)*Ex; 
  var = mu * (1 + mu / exp(log_dispersion));
  nll -= dnbinom2(dx, mu, var, 1).sum();

  matrix<Type> LQ_params_mat(4, n_periods);
  for(int i=0; i<n_periods; i++){
  LQ_params_mat(0,i) = Type(1.0);
  }
  LQ_params_mat.row(1) = h_params;
  LQ_params_mat.row(2) = h2_params;
  LQ_params_mat.row(3) = k_params;

  matrix<Type> log_mx_mat = LQ_baseline*LQ_params_mat;
  matrix<Type> mx_mat = exp(log_mx_mat.array());
  matrix<Type> weights_mat(mx_mat.rows() - open_idx + 1, n_periods);
  weights_mat = 1 / (1 + 2.5*mx_mat.block(open_idx - 1, 0, mx_mat.rows() - open_idx + 1, n_periods).array());
  
  for(int i=0; i < n_periods; i++) {
  	for(int j=open_idx - 1; j < mx_mat.rows()-1; j++) {
		for(int k=j-open_idx+2; k < mx_mat.rows()-open_idx+1; k++) {
  			weights_mat(k,i) *= (1 - 2.5*mx_mat(j,i)) / (1 + 2.5*mx_mat(j,i));
		}
  	}
  }

  vector<Type> weights_sum = weights_mat.colwise().sum();
  
  for(int i=0; i < n_periods; i++) {
    for(int j=0; j < weights_mat.rows(); j++) {
    weights_mat(j,i) *= 1 / weights_sum(i);
    }
  }
  
  matrix<Type> mx_open_weighted = mx_mat.block(open_idx - 1, 0, mx_mat.rows() - open_idx + 1, n_periods).array() * weights_mat.array();

  matrix<Type> mx_mat_aggr(open_idx, n_periods);
  mx_mat_aggr.block(0, 0, open_idx - 1, n_periods) = mx_mat.block(0, 0, open_idx - 1, n_periods);
  mx_mat_aggr.row(open_idx - 1) = mx_open_weighted.colwise().sum();
  
  matrix<Type> sx_mat(mx_mat_aggr.rows() + 1, n_periods);
  
  for(int i=0; i < n_periods; i++) {
    //UDD 0.5q0 (2.5q[0-5] in this case)
    sx_mat(0,i) = 1 / (1 + 2.5 * mx_mat_aggr(0,i)); 
    //UDD open age group sx
    sx_mat(mx_mat_aggr.rows(), i) = (1 - 2.5*mx_mat_aggr(mx_mat_aggr.rows() - 1, i)) / ((1 + 2.5*mx_mat_aggr(mx_mat_aggr.rows() - 1, i)));
    
    for(int j=0; j < (mx_mat_aggr.rows() - 1); j++) {
      sx_mat(j+1,i) = (1 - 2.5*mx_mat_aggr(j,i)) / (1 + 2.5*mx_mat_aggr(j+1,i));
     }
   }

  // prior for log(fx)
  nll -= dnorm(log_fx, log_fx_mean, sigma_fx, true).sum();
  vector<Type> fx(exp(log_fx));
  MapMatrixXXt fx_mat(fx.data(), n_fx, n_periods);

  // prior for gx
  nll -= dnorm(gx, gx_mean, sigma_gx, true).sum();
  MapMatrixXXt gx_mat(gx.data(), basepop.size(), n_periods);

  // population projection
  PopulationProjection<Type> proj(ccmpp<Type>(basepop, sx_mat, fx_mat, gx_mat,
					      srb, interval, fx_idx-1));
  
  // likelihood for log census counts
  for(int i = 0; i < census_year_idx.size(); i++) {
    nll -= dnorm(vector<Type>(census_log_pop.col(i)),
  		 log(vector<Type>(proj.population.col(census_year_idx[i] - 1))),
  		 sigma_logpop, true).sum();
  }


  MatrixXXt proj_period_deaths(proj.period_deaths());
  vector<Type> period_deaths(MapVectorXt(proj_period_deaths.data(),
					 proj_period_deaths.size()));
  vector<Type> births(MapVectorXt(proj.births.data(),
  				  proj.births.size()));

  DATA_INTEGER(calc_outputs);

  if(calc_outputs) {
    
    vector<Type> population(MapVectorXt(proj.population.data(),
					proj.population.size()));
    vector<Type> cohort_deaths(MapVectorXt(proj.cohort_deaths.data(),
					   proj.cohort_deaths.size()));
    vector<Type> infants(MapVectorXt(proj.infants.data(),
				     proj.infants.size()));
    vector<Type> migrations(MapVectorXt(proj.migrations.data(),
					proj.migrations.size()));
    REPORT(population);
    REPORT(cohort_deaths);
    REPORT(period_deaths);
    REPORT(births);
    REPORT(infants);
    REPORT(migrations);
    REPORT(sx_mat);
    REPORT(mx_mat);
    REPORT(fx);
    REPORT(gx);
  }
  
  return Type(nll);

}

