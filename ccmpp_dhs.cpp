#ifndef CCMPP_VR_TMB_HPP
#define CCMPP_VR_TMB_HPP

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

#include "../ccmpp.h"

template<class Type>
Type ccmpp_vr_tmb(objective_function<Type>* obj)
{
  std::cout << "ccmpp_vr_tmb" << std::endl;
  
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
  DATA_SCALAR(interval);
  DATA_INTEGER(n_periods);
  DATA_INTEGER(fx_idx);
  DATA_INTEGER(n_fx);
  //DATA_VECTOR(log_hiv_mort_mean);
  //DATA_INTEGER(hiv_mort_start);
  //DATA_INTEGER(hiv_mort_end);

  // census data
  DATA_MATRIX(census_log_pop);
  DATA_IVECTOR(census_year_idx);

  Type nll(0.0);

  // Hyper priors
  PARAMETER(log_tau2_logpop);
  nll -= dlgamma(log_tau2_logpop, Type(1.0), Type(1.0 / 0.0109), true);
  Type sigma_logpop(exp(-0.5 * log_tau2_logpop));

  PARAMETER(log_tau2_sx);
  nll -= dlgamma(log_tau2_sx, Type(1.0), Type(1.0 / 0.0109), true);
  Type sigma_sx(exp(-0.5 * log_tau2_sx));

  PARAMETER(log_tau2_fx);
  nll -= dlgamma(log_tau2_fx, Type(1.0), Type(1.0 / 0.0109), true);
  Type sigma_fx(exp(-0.5 * log_tau2_fx));

  PARAMETER(log_tau2_gx);
  nll -= dlgamma(log_tau2_gx, Type(1.0), Type(1.0 / 0.0436), true);
  Type sigma_gx(exp(-0.5 * log_tau2_gx));
  
  PARAMETER(log_lambda_tp);
  nll -= dnorm(log_lambda_tp,Type(0.0),Type(5.0),1);
  
  // prior for base population
  PARAMETER_VECTOR(log_basepop);
  nll -= dnorm(log_basepop, log_basepop_mean, sigma_logpop, true).sum();
  vector<Type> basepop(exp(log_basepop));
  
  // prior for log(mx)
  DATA_VECTOR(mx_params_mean);
  DATA_VECTOR(log_baseline_mx);
  DATA_VECTOR(baseline_mx_DX);
  DATA_SPARSE_MATRIX(mx_params_DX);
  DATA_SPARSE_MATRIX(tp_DX);
  DATA_SPARSE_MATRIX(penal_tp);
  DATA_SPARSE_MATRIX(null_penal_tp);
  PARAMETER_VECTOR(mx_params);
  PARAMETER_VECTOR(tp_params);
  PARAMETER(log_dispersion);

  nll -= dnorm(mx_params, mx_params_mean, sigma_sx, true).sum();
  nll -= dnorm(log_dispersion,Type(0.0),Type(5.0),1);
  SparseMatrix<Type> QQ_tp = exp(log_lambda_tp)*penal_tp+null_penal_tp;
  nll += GMRF(QQ_tp)(tp_params);

  //likelihood for DHS data
  vector<Type> mu, var;
  DATA_VECTOR(dx);
  DATA_VECTOR(Ex);
  mu = exp(baseline_mx_DX + mx_params_DX*mx_params + tp_DX*tp_params)*Ex; 
  var = mu * (1 + mu / exp(log_dispersion));

  // likelihood for deaths
  nll -= dnbinom2(dx, mu, var, 1).sum();

 // PARAMETER(log_tau2_hiv_mort);
 // PARAMETER_VECTOR(log_hiv_mort);
 // CHANGE PRIOR FOR HIV MORT, MIGHT BE TOO FLAT AT LOW HIV MORT RANGE

 // nll -= dlgamma(log_tau2_hiv_mort, Type(1.0), Type(1.0 / 0.0436), true);
 // Type sigma_hiv_mort(exp(-0.5 * log_tau2_hiv_mort));
 // nll -= dnorm(log_hiv_mort, log_hiv_mort_mean, sigma_hiv_mort, true).sum();   

 
  // matrix<Type> hiv_mort_mat(basepop.size()+1, n_periods);
  
  // for(int i = 0; i < n_periods; i++) {
  // for(int j = 0; j < basepop.size()+1; j++) {
  //   hiv_mort_mat(j,i) = Type(1.0);
  //    }
  // }
  
  //for(int i = 0; i < n_periods; i++) {
  //  for(int j = hiv_mort_start-1; j < hiv_mort_end; j++) {
  //    hiv_mort_mat(j,i) = exp(-exp(log_hiv_mort(i)));
  //    }
  // }

  //create sx_mat for ccmpp
  vector<Type> log_mx(log_baseline_mx.size()*n_periods);
  for(int i = 0; i < n_periods; i++) {
    log_mx.segment(i*log_baseline_mx.size(),log_baseline_mx.size()) = log_baseline_mx + mx_params(i);
  }
  vector<Type> mx(exp(log_mx));
  MapMatrixXXt mx_mat(mx.data(), basepop.size(), n_periods);
  
  matrix<Type> sx_mat(basepop.size() + 1, n_periods);
  
  for(int i=0; i < n_periods; i++) {
    //UDD 0.5q0 (2.5q[0-5] in this case)
    sx_mat(0,i) = 1 / (1 + 0.5 * mx_mat(0,i)); 
    //UDD open age group sx
    sx_mat(log_baseline_mx.size(), i) = (1 - 0.5*mx_mat(log_baseline_mx.size() - 1, i)) / ((1 + 0.5*mx_mat(log_baseline_mx.size() - 1, i)));
    
    for(int j=0; j < (log_baseline_mx.size() - 1); j++) {
      sx_mat(j+1,i) = (1 - 0.5*mx_mat(j,i)) / (1 + 0.5*mx_mat(j+1,i));
    }
  }
  
  
 // matrix<Type> all_mort_mat(basepop.size()+1, n_periods);
 // all_mort_mat = sx_mat.array() * hiv_mort_mat.array();
  
  // prior for log(fx)
  PARAMETER_VECTOR(log_fx);
  nll -= dnorm(log_fx, log_fx_mean, sigma_fx, true).sum();
  vector<Type> fx(exp(log_fx));
  MapMatrixXXt fx_mat(fx.data(), n_fx, n_periods);

  // prior for gx
  PARAMETER_VECTOR(gx);
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
    REPORT(fx);
    REPORT(gx);
  }
  
  return Type(nll);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif