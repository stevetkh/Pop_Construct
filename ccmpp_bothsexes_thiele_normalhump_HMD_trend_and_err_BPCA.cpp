#include <TMB.hpp>                                
#include "ccmpp.h"
#include "ccmpp_m.h"

template <typename Type>
vector<Type> sumByAG1D(const vector<Type>& B, const vector<int>& pace_v, int new_row)
{
  vector<Type> A(new_row);
  A.setZero();
  int r_i = pace_v(0); // first age group
  for (int i = 0; i < B.size(); ++i) {
    if (pace_v(i) != r_i)
      ++r_i;
    A(r_i - 1) += B(i);
  } // end age-groups
  return A;
}

template <typename Type>
Type gumbel_density(const Type& tau, const Type& theta)
{
  Type ll = 0.0;
  ll += log(theta) - log(Type(2.0)) - Type(1.5) * log(tau) - theta / sqrt(tau);
  return ll;
}

template <typename Type>
Type get_first_var_2d(const vector<Type>& A, const vector<Type>& cov1, const vector<Type>& cov2)
{
  int s = A.size();
  Type v = 0.0;

  for(int i = 0; i < s; i++){
    for(int j = 0; j < s; j++){
      for(int k = 0; k < s; k++){
        for(int l = 0; l < s; l++){
          v += A(i) * A(j) * A(k) * A(l) * cov1(abs(i-j)) * cov2(abs(k-l));    
        }
      }
    }
  } 
  return v;
}

template <typename Type>
Type AR2_rho_prior_base1(const Type& rho, const Type& theta)
{
  Type ll = 0.0;
  ll += -theta * sqrt((Type(1.0) - rho) * (Type(1.0) - rho) / (Type(1.0) + rho)) + log(Type(3.0) + rho) - Type(1.5) * log(Type(1.0) + rho);
// ll += log(theta) - log(Type(2.0)) - log(Type(1.0) - exp(-theta));
  return ll;
}

template <typename Type>
Type AR1_rho_prior_base0(const Type& rho, const Type& theta)
{
  Type ll = 0.0;
  ll += log(theta) - log(Type(2.0)) - theta * sqrt(-log(1.0 - rho * rho)) + log(rho) - log(1 - rho * rho) - 0.5 * log(-log(1 - rho * rho));
  return ll;
}

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
  DATA_INTEGER(n_periods);
  DATA_INTEGER(fx_idx);
  DATA_INTEGER(n_fx);
  DATA_IVECTOR(pop_start);
  DATA_IVECTOR(pop_end);
  DATA_INTEGER(open_idx);

  DATA_MATRIX(census_log_pop_f_5);
  DATA_MATRIX(census_log_pop_m_5);
  DATA_IVECTOR(census_year_idx_5);
  DATA_IVECTOR(census_year_group_idx_5);
  DATA_INTEGER(n_agegrp_5);
  DATA_IVECTOR(pop_start_5);
  DATA_IVECTOR(pop_end_5);
  DATA_INTEGER(calc_5q0);

  PARAMETER_VECTOR(log_tau2_logpop); //1 = WPP males, 2 = data males, 3 = WPP females, 4 = data females
  PARAMETER(log_marginal_lambda_fx);
  PARAMETER(log_marginal_lambda_gx); //1 = male, 2 = female

  PARAMETER_VECTOR(log_basepop_f);
  PARAMETER_VECTOR(log_basepop_m);
  PARAMETER_VECTOR(log_fx_spline_params);
  PARAMETER_VECTOR(gx_f_spline_params);
  PARAMETER_VECTOR(gx_m_spline_params);

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

  DATA_VECTOR(df_br);
  DATA_VECTOR(dm_br);
  DATA_VECTOR(Ef_br);
  DATA_VECTOR(Em_br);
  DATA_IVECTOR(df_age_br);
  DATA_IVECTOR(dm_age_br);
  DATA_IVECTOR(df_time_br);
  DATA_IVECTOR(dm_time_br);

  DATA_VECTOR(thiele_age);
  
  DATA_SPARSE_MATRIX(D_time);
  DATA_SPARSE_MATRIX(D_err);
  DATA_VECTOR(slope_ind);
  
  DATA_SPARSE_MATRIX(tp_mat_br_f_age0);
  DATA_SPARSE_MATRIX(tp_mat_br_m_age0);
  DATA_SPARSE_MATRIX(tp_mat_br_f_age1);
  DATA_SPARSE_MATRIX(tp_mat_br_m_age1);
  DATA_SPARSE_MATRIX(tp_mat_br_f_age5);
  DATA_SPARSE_MATRIX(tp_mat_br_m_age5);
  DATA_SPARSE_MATRIX(tp_mat_br_f_age10);
  DATA_SPARSE_MATRIX(tp_mat_br_m_age10);

  DATA_SPARSE_MATRIX(penal_tp);
  DATA_SPARSE_MATRIX(penal_tp_0);
  DATA_SPARSE_MATRIX(null_penal_tp);
  DATA_VECTOR(tp_mean);
  DATA_SPARSE_MATRIX(null_penal_tp_br);
  DATA_VECTOR(tp_mean_br);
  
  PARAMETER(log_dispersion); //1 = male, 2 = female
  PARAMETER(log_dispersion_br); //1 = male, 2 = female

  PARAMETER(log_lambda_tp);
  PARAMETER_VECTOR(tp_params);
  PARAMETER(tp_slope);
  PARAMETER(tp_params_5);
  PARAMETER(tp_params_10);

  PARAMETER(log_lambda_tp_br_age0);
  PARAMETER_VECTOR(tp_params_br_age0);
  PARAMETER(tp_slope_br_age0);
  PARAMETER(tp_params_5_br_age0);
  PARAMETER(tp_params_10_br_age0);

  PARAMETER(log_lambda_tp_br_age1);
  PARAMETER_VECTOR(tp_params_br_age1);
  PARAMETER(tp_slope_br_age1);
  PARAMETER(tp_params_5_br_age1);
  PARAMETER(tp_params_10_br_age1);

  PARAMETER(log_lambda_tp_br_age5);
  PARAMETER_VECTOR(tp_params_br_age5);
  PARAMETER(tp_slope_br_age5);
  PARAMETER(tp_params_5_br_age5);
  PARAMETER(tp_params_10_br_age5);

  PARAMETER(log_lambda_tp_br_age10);
  PARAMETER_VECTOR(tp_params_br_age10);
  PARAMETER(tp_slope_br_age10);
  PARAMETER(tp_params_5_br_age10);
  PARAMETER(tp_params_10_br_age10);

  DATA_SCALAR(theta_tp);

  DATA_VECTOR(log_spline_prior_means_m);  
  DATA_VECTOR(log_spline_prior_means_f);  

  PARAMETER_VECTOR(log_spline_means_m);
  PARAMETER_VECTOR(log_spline_means_f);
  PARAMETER_VECTOR(log_spline_slopes_m);
  PARAMETER_VECTOR(log_spline_slopes_f);
  PARAMETER_MATRIX(log_spline_err_m);
  PARAMETER_MATRIX(log_spline_err_f);

  PARAMETER_MATRIX(log_spline_hump_m);
  PARAMETER_MATRIX(log_spline_hump_f);

  PARAMETER_VECTOR(log_prec_all_RW2_hump);

  DATA_INTEGER(loghump);

  PARAMETER(logit_rho_delta);
  PARAMETER(logit_rho_epsilon);
  DATA_SCALAR(theta_rho_hump);

  DATA_MATRIX(V_mf_mean);
  DATA_MATRIX(V_mf_slope);
  DATA_MATRIX(V_mf);
  
  DATA_MATRIX(PCA);  
  
  PARAMETER(log_V_factor_mean);
  PARAMETER(log_V_factor_slope);
  PARAMETER(log_V_factor);
  
  DATA_VECTOR(log_V_factor_prior_sd);
  
  DATA_SCALAR(lambda_prior_sd);
  DATA_SCALAR(delta_prior_sd);
  DATA_SCALAR(epsilon_prior_sd);
  
  DATA_SPARSE_MATRIX(D_agetime);
  DATA_SPARSE_MATRIX(D_agetime_fert);
  DATA_VECTOR(D_firstrow);  
  DATA_VECTOR(D_firstrow_ARIMA);

  PARAMETER(logit_rho_fx_age);
  PARAMETER(logit_rho_fx_time);
  PARAMETER(logit_rho_gx_age);
  PARAMETER(logit_rho_gx_time);

  DATA_SCALAR(upper_marginal_sd_fx);
  DATA_SCALAR(upper_marginal_sd_gx);

  // DATA_SCALAR(theta_logpop);

  DATA_SCALAR(theta_rho_fx_age_base1);
  DATA_SCALAR(theta_rho_fx_time_base1);
  DATA_SCALAR(theta_rho_gx_age_base1);
  DATA_SCALAR(theta_rho_gx_time_base1);

  DATA_MATRIX(census_deaths_f);
  DATA_MATRIX(census_deaths_m);
  DATA_IVECTOR(census_deaths_year_idx);
  DATA_IVECTOR(deaths_start);
  DATA_IVECTOR(deaths_end);

  DATA_MATRIX(census_deaths_f_5);
  DATA_MATRIX(census_deaths_m_5);
  DATA_IVECTOR(census_deaths_year_idx_5);
  DATA_IVECTOR(census_deaths_group_idx_5);
  DATA_INTEGER(n_agegrp_5_deaths);
  DATA_IVECTOR(deaths_start_5);
  DATA_IVECTOR(deaths_end_5);

  PARAMETER_VECTOR(log_dispersion_census_deaths);
  PARAMETER_VECTOR(log_dispersion_census_deaths_5);

  PARAMETER(log_prec_5q0);
  DATA_SCALAR(theta_5q0);

  DATA_VECTOR(log_dat_v5q0_m);
  DATA_VECTOR(log_dat_v5q0_f);
  DATA_IVECTOR(year_ind_5q0);

  Type nll(0.0);
  
  int C_dim = log_spline_err_m.cols();
  int no_basis_time_int = log_spline_hump_m.rows();

  Type lambda_tp = exp(log_lambda_tp);
  nll -= gumbel_density(lambda_tp, theta_tp);
  nll -= log_lambda_tp;

  Type lambda_tp_br_age0 = exp(log_lambda_tp_br_age0);
  nll -= gumbel_density(lambda_tp_br_age0, theta_tp);
  nll -= log_lambda_tp_br_age0;

  Type lambda_tp_br_age1 = exp(log_lambda_tp_br_age1);
  nll -= gumbel_density(lambda_tp_br_age1, theta_tp);
  nll -= log_lambda_tp_br_age1;

  Type lambda_tp_br_age5 = exp(log_lambda_tp_br_age5);
  nll -= gumbel_density(lambda_tp_br_age5, theta_tp);
  nll -= log_lambda_tp_br_age5;

  Type lambda_tp_br_age10 = exp(log_lambda_tp_br_age10);
  nll -= gumbel_density(lambda_tp_br_age10, theta_tp);
  nll -= log_lambda_tp_br_age10;

  SparseMatrix<Type> QQ_tp = lambda_tp * penal_tp + 4 * penal_tp_0 + null_penal_tp;
  nll += GMRF(QQ_tp)(tp_params - tp_slope * tp_mean);
  nll -= dnorm(tp_slope, Type(0.0), Type(0.1), true);
  nll -= dnorm(tp_params_5, Type(0.0), Type(0.1), true);
  nll -= dnorm(tp_params_10, Type(0.0), Type(0.1), true);
  tp_params(5) += tp_params_5;
  tp_params(10) += tp_params_10;

  SparseMatrix<Type> QQ_tp_br_age0 = lambda_tp_br_age0 * penal_tp + 4 * penal_tp_0 + null_penal_tp_br;
  nll += GMRF(QQ_tp_br_age0)(tp_params_br_age0 - tp_slope_br_age0 * tp_mean_br);
  nll -= dnorm(tp_slope_br_age0, Type(0.0), Type(0.1), true);
  nll -= dnorm(tp_params_5_br_age0, Type(0.0), Type(0.1), true);
  nll -= dnorm(tp_params_10_br_age0, Type(0.0), Type(0.1), true);
  tp_params_br_age0(5) += tp_params_5_br_age0;
  tp_params_br_age0(10) += tp_params_10_br_age0;

  SparseMatrix<Type> QQ_tp_br_age1 = lambda_tp_br_age1 * penal_tp + 4 * penal_tp_0 + null_penal_tp_br;
  nll += GMRF(QQ_tp_br_age1)(tp_params_br_age1 - tp_slope_br_age1 * tp_mean_br);
  nll -= dnorm(tp_slope_br_age1, Type(0.0), Type(0.1), true);
  nll -= dnorm(tp_params_5_br_age1, Type(0.0), Type(0.1), true);
  nll -= dnorm(tp_params_10_br_age1, Type(0.0), Type(0.1), true);
  tp_params_br_age1(5) += tp_params_5_br_age1;
  tp_params_br_age1(10) += tp_params_10_br_age1;

  SparseMatrix<Type> QQ_tp_br_age5 = lambda_tp_br_age5 * penal_tp + 4 * penal_tp_0 + null_penal_tp_br;
  nll += GMRF(QQ_tp_br_age5)(tp_params_br_age5 - tp_slope_br_age5 * tp_mean_br);
  nll -= dnorm(tp_slope_br_age5, Type(0.0), Type(0.1), true);
  nll -= dnorm(tp_params_5_br_age5, Type(0.0), Type(0.1), true);
  nll -= dnorm(tp_params_10_br_age5, Type(0.0), Type(0.1), true);
  tp_params_br_age5(5) += tp_params_5_br_age5;
  tp_params_br_age5(10) += tp_params_10_br_age5;

  SparseMatrix<Type> QQ_tp_br_age10 = lambda_tp_br_age10 * penal_tp + 4 * penal_tp_0 + null_penal_tp_br;
  nll += GMRF(QQ_tp_br_age10)(tp_params_br_age10 - tp_slope_br_age10 * tp_mean_br);
  nll -= dnorm(tp_slope_br_age10, Type(0.0), Type(0.1), true);
  nll -= dnorm(tp_params_5_br_age10, Type(0.0), Type(0.1), true);
  nll -= dnorm(tp_params_10_br_age10, Type(0.0), Type(0.1), true);
  tp_params_br_age10(5) += tp_params_5_br_age10;
  tp_params_br_age10(10) += tp_params_10_br_age10;

  vector<Type> rho_vec(2);
  rho_vec(0) = 2.0;
  rho_vec(1) = -1.0;

  vector<Type> sigma_RW2_hump = exp(-0.5 * log_prec_all_RW2_hump);
  nll -= dnorm(sigma_RW2_hump(0), Type(0), lambda_prior_sd, true);
  nll -= dnorm(sigma_RW2_hump(1), Type(0), delta_prior_sd, true);
  nll -= dnorm(sigma_RW2_hump(2), Type(0), epsilon_prior_sd, true);
  nll -= -Type(0.5) * log_prec_all_RW2_hump.sum();
  
  Type rho_delta = invlogit(logit_rho_delta);
  nll -= AR2_rho_prior_base1(rho_delta, theta_rho_hump);
  nll -= Type(2.0) * log(rho_delta) - logit_rho_delta;

  Type rho_epsilon = invlogit(logit_rho_epsilon);
  nll -= AR2_rho_prior_base1(rho_epsilon, theta_rho_hump);
  nll -= Type(2.0) * log(rho_epsilon) - logit_rho_epsilon;

  nll -= dnorm(log_V_factor_mean, Type(0.0), Type(log_V_factor_prior_sd(0)), true);
  nll -= dnorm(log_V_factor_slope, Type(0.0), Type(log_V_factor_prior_sd(1)), true);  
  nll -= dnorm(log_V_factor, Type(0.0), Type(log_V_factor_prior_sd(2)), true);

  matrix<Type> V_mf_mean_adjusted(V_mf_mean.rows(), V_mf_mean.cols());
  matrix<Type> V_mf_slope_adjusted(V_mf_slope.rows(), V_mf_slope.cols());
  matrix<Type> V_mf_adjusted(V_mf.rows(), V_mf.cols());
  V_mf_mean_adjusted = exp(log_V_factor_mean) * V_mf_mean;
  V_mf_slope_adjusted = exp(log_V_factor_slope) * V_mf_slope;
  V_mf_adjusted = exp(log_V_factor) * V_mf;
  
  MVNORM_t<Type> MVN_mean_mf(V_mf_mean_adjusted);
  MVNORM_t<Type> MVN_slope_mf(V_mf_slope_adjusted);
  MVNORM_t<Type> MVN_err_mf(V_mf_adjusted);
    
  vector<Type> log_spline_means_mf_demean_nohump(C_dim * 2);
  log_spline_means_mf_demean_nohump(0) = log_spline_means_m(0) - log_spline_prior_means_m(0);
  log_spline_means_mf_demean_nohump(1) = log_spline_means_m(1) - log_spline_prior_means_m(1);
  log_spline_means_mf_demean_nohump(2) = log_spline_means_m(2) - log_spline_prior_means_m(2);
  log_spline_means_mf_demean_nohump(3) = log_spline_means_m(3) - log_spline_prior_means_m(6);
  log_spline_means_mf_demean_nohump(4) = log_spline_means_m(4) - log_spline_prior_means_m(7);  
  log_spline_means_mf_demean_nohump(5) = log_spline_means_f(0) - log_spline_prior_means_f(0);
  log_spline_means_mf_demean_nohump(6) = log_spline_means_f(1) - log_spline_prior_means_f(1);
  log_spline_means_mf_demean_nohump(7) = log_spline_means_f(2) - log_spline_prior_means_f(2);
  log_spline_means_mf_demean_nohump(8) = log_spline_means_f(3) - log_spline_prior_means_f(6);
  log_spline_means_mf_demean_nohump(9) = log_spline_means_f(4) - log_spline_prior_means_f(7);  
  nll += MVN_mean_mf(log_spline_means_mf_demean_nohump);

  vector<Type> temp_vec_slopes(C_dim * 2);
  temp_vec_slopes << log_spline_slopes_m, log_spline_slopes_f;
  nll += MVN_slope_mf(temp_vec_slopes);
  
  vector<Type> temp_vec_err_m(C_dim);
  vector<Type> temp_vec_err_f(C_dim);
  vector<Type> temp_vec_err(C_dim * 2);
  for(int t = 0; t < D_err.cols(); t++){
	  temp_vec_err_m = log_spline_err_m.row(t);
	  temp_vec_err_f = log_spline_err_f.row(t);
	  temp_vec_err << temp_vec_err_m, temp_vec_err_f;
	  nll += MVN_err_mf(temp_vec_err);
    }

  nll -= dnorm(diff(diff(vector<Type>(log_spline_hump_m.col(0)))), Type(0.0), Type(sigma_RW2_hump(0)), true).sum();
  nll -= dnorm(diff(diff(vector<Type>(log_spline_hump_f.col(0)))), Type(0.0), Type(sigma_RW2_hump(0)), true).sum();
  nll += SCALE(ARk(vector<Type>(rho_delta * rho_vec)), sigma_RW2_hump(1))(vector<Type>(log_spline_hump_m.col(1)) - log_spline_prior_means_m(4));
  nll += SCALE(ARk(vector<Type>(rho_delta * rho_vec)), sigma_RW2_hump(1))(vector<Type>(log_spline_hump_f.col(1)) - log_spline_prior_means_f(4));
  nll += SCALE(ARk(vector<Type>(rho_epsilon * rho_vec)), sigma_RW2_hump(2))(vector<Type>(log_spline_hump_m.col(2)) - log_spline_prior_means_m(5));
  nll += SCALE(ARk(vector<Type>(rho_epsilon * rho_vec)), sigma_RW2_hump(2))(vector<Type>(log_spline_hump_f.col(2)) - log_spline_prior_means_f(5));
  
  array<Type> log_par_m_temp(D_err.rows(), C_dim);
  array<Type> log_par_f_temp(D_err.rows(), C_dim); 
  array<Type> log_par_m_hump(D_err.rows(), 3);
  array<Type> log_par_f_hump(D_err.rows(), 3); 
  
  log_par_m_temp = D_err * log_spline_err_m;
  log_par_f_temp = D_err * log_spline_err_f;
  log_par_m_hump = D_time * log_spline_hump_m;
  log_par_f_hump = D_time * log_spline_hump_f;
	
	for(int p = 0; p < C_dim; p++){
		log_par_m_temp.col(p) += log_spline_means_m(p) + log_spline_slopes_m(p) * slope_ind;
		log_par_f_temp.col(p) += log_spline_means_f(p) + log_spline_slopes_f(p) * slope_ind;
		}
  
  array<Type> log_par_m(D_err.rows(), 8);
  array<Type> log_par_f(D_err.rows(), 8);
  log_par_m.col(0) = log_par_m_temp.col(0);
  log_par_m.col(1) = log_par_m_temp.col(1);
  log_par_m.col(2) = log_par_m_temp.col(2);
  log_par_m.col(3) = log_par_m_hump.col(0);
  log_par_m.col(4) = log_par_m_hump.col(1);
  log_par_m.col(5) = log_par_m_hump.col(2);
  log_par_m.col(6) = log_par_m_temp.col(3);
  log_par_m.col(7) = log_par_m_temp.col(4);
  
  log_par_f.col(0) = log_par_f_temp.col(0);
  log_par_f.col(1) = log_par_f_temp.col(1);
  log_par_f.col(2) = log_par_f_temp.col(2);
  log_par_f.col(3) = log_par_f_hump.col(0);
  log_par_f.col(4) = log_par_f_hump.col(1);
  log_par_f.col(5) = log_par_f_hump.col(2);
  log_par_f.col(6) = log_par_f_temp.col(3);
  log_par_f.col(7) = log_par_f_temp.col(4);
  
  vector<Type> phi_f_vec = log_par_f.col(0);
  vector<Type> disp_f_vec = log_par_f.col(1);
  vector<Type> psi_f_vec = log_par_f.col(2);
  vector<Type> lambda_f_vec = log_par_f.col(3);
  vector<Type> delta_f_vec = log_par_f.col(4);
  vector<Type> epsilon_f_vec = log_par_f.col(5);
  vector<Type> A_f_vec = log_par_f.col(6);
  vector<Type> B_f_vec = log_par_f.col(7);

  vector<Type> phi_m_vec = log_par_m.col(0);
  vector<Type> disp_m_vec = log_par_m.col(1);
  vector<Type> psi_m_vec = log_par_m.col(2);
  vector<Type> lambda_m_vec = log_par_m.col(3);
  vector<Type> delta_m_vec = log_par_m.col(4);
  vector<Type> epsilon_m_vec = log_par_m.col(5);
  vector<Type> A_m_vec = log_par_m.col(6);
  vector<Type> B_m_vec = log_par_m.col(7);
 
  matrix<Type> phi_f = exp(phi_f_vec);
  matrix<Type> disp_f = exp(disp_f_vec);
  matrix<Type> psi_f = exp(psi_f_vec);
  matrix<Type> lambda_f = exp(lambda_f_vec);
  matrix<Type> delta_f = exp(delta_f_vec);
  matrix<Type> epsilon_f = exp(epsilon_f_vec);
  matrix<Type> A_f = A_f_vec;
  matrix<Type> B_f = B_f_vec;

  matrix<Type> phi_m = exp(phi_m_vec);
  matrix<Type> disp_m = exp(disp_m_vec);
  matrix<Type> psi_m = exp(psi_m_vec);
  matrix<Type> lambda_m = exp(lambda_m_vec);
  matrix<Type> delta_m = exp(delta_m_vec);
  matrix<Type> epsilon_m = exp(epsilon_m_vec);
  matrix<Type> A_m = A_m_vec;
  matrix<Type> B_m = B_m_vec;

  matrix<Type> mx_mat_f(thiele_age.size(), D_time.rows()); 
  matrix<Type> mx_mat_m(thiele_age.size(), D_time.rows()); 
  vector<Type> log_v5q0_m(phi_m.size());
  vector<Type> log_v5q0_f(phi_f.size());
  vector<Type> PCA_one = PCA.col(0);
  vector<Type> PCA_two = PCA.col(1);
  
  if(loghump) {
      for(int i = 0; i < D_err.rows(); i++){
		  mx_mat_f.col(i) = pow(phi_f(i), vector<Type> (pow(thiele_age + disp_f(i), psi_f(i)))) + 
			lambda_f(i)*exp(-delta_f(i)*((log(thiele_age) - log(epsilon_f(i)))*(log(thiele_age) - log(epsilon_f(i))))) +
			exp(A_f(i) * PCA_one + B_f(i) * PCA_two);
		  mx_mat_m.col(i) = pow(phi_m(i), vector<Type> (pow(thiele_age + disp_m(i), psi_m(i)))) + 
			lambda_m(i)*exp(-delta_m(i)*((log(thiele_age) - log(epsilon_m(i)))*(log(thiele_age) - log(epsilon_m(i))))) +
			exp(A_m(i) * PCA_one + B_m(i) * PCA_two);
						
		  log_v5q0_m(i) = log(Type(1.0) - vector<Type>(Type(1.0) - vector<Type>(mx_mat_m.col(i).head(5)) / (Type(1.0) + Type(0.5) * vector<Type>(mx_mat_m.col(i).head(5)))).prod());
		  log_v5q0_f(i) = log(Type(1.0) - vector<Type>(Type(1.0) - vector<Type>(mx_mat_f.col(i).head(5)) / (Type(1.0) + Type(0.5) * vector<Type>(mx_mat_f.col(i).head(5)))).prod());  
      }
  } else {
      for(int i = 0; i < D_err.rows(); i++){
		  mx_mat_f.col(i) = pow(phi_f(i), vector<Type> (pow(thiele_age + disp_f(i), psi_f(i)))) + 
			lambda_f(i)*exp(-delta_f(i)*((thiele_age - epsilon_f(i))*(thiele_age - epsilon_f(i)))) + 
			exp(A_f(i) * PCA_one + B_f(i) * PCA_two);
		  mx_mat_m.col(i) = pow(phi_m(i), vector<Type> (pow(thiele_age + disp_m(i), psi_m(i)))) + 
			lambda_m(i)*exp(-delta_m(i)*((thiele_age - epsilon_m(i))*(thiele_age - epsilon_m(i)))) +
			exp(A_m(i) * PCA_one + B_m(i) * PCA_two);

		  log_v5q0_m(i) = log(Type(1.0) - vector<Type>(Type(1.0) - vector<Type>(mx_mat_m.col(i).head(5)) / (Type(1.0) + Type(0.5) * vector<Type>(mx_mat_m.col(i).head(5)))).prod());
		  log_v5q0_f(i) = log(Type(1.0) - vector<Type>(Type(1.0) - vector<Type>(mx_mat_f.col(i).head(5)) / (Type(1.0) + Type(0.5) * vector<Type>(mx_mat_f.col(i).head(5)))).prod());  
      } 
  }
 
  //fit to 5q0
  Type lambda_5q0 = exp(log_prec_5q0);
  nll -= gumbel_density(lambda_5q0, theta_5q0);
  nll -= log_prec_5q0;

  Type sigma_5q0 = exp(-0.5 * log_prec_5q0);

  if(calc_5q0){
	  for(int i = 0; i < year_ind_5q0.size(); i++){
		  nll -= dnorm(log_dat_v5q0_m(i), log_v5q0_m(year_ind_5q0(i)), sigma_5q0, 1) + dnorm(log_dat_v5q0_f(i), log_v5q0_f(year_ind_5q0(i)), sigma_5q0, 1);
	  }
	}

  //likelihood for DHS data
  nll -= dnorm(log_dispersion, Type(0.0), Type(10.0), 1);

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

  varf = muf * (1 + muf / exp(log_dispersion));
  varm = mum * (1 + mum / exp(log_dispersion));
  nll -= dnbinom2(df, muf, varf, 1).sum();
  nll -= dnbinom2(dm, mum, varm, 1).sum();

  //likelihood for BR data
  nll -= dnorm(log_dispersion_br, Type(0.0), Type(10.0), 1);

  vector<Type> mum_br_tp(dm_br.size());
  vector<Type> muf_br_tp(df_br.size());
 
  mum_br_tp = tp_mat_br_m_age0 * tp_params_br_age0 + tp_mat_br_m_age1 * tp_params_br_age1 +
	      tp_mat_br_m_age5 * tp_params_br_age5 + tp_mat_br_m_age10 * tp_params_br_age10;
  
  muf_br_tp = tp_mat_br_f_age0 * tp_params_br_age0 + tp_mat_br_f_age1 * tp_params_br_age1 +
	      tp_mat_br_f_age5 * tp_params_br_age5 + tp_mat_br_f_age10 * tp_params_br_age10;

  vector<Type> muf_br(df_br.size());
  vector<Type> varf_br(df_br.size());
  vector<Type> mum_br(dm_br.size());
  vector<Type> varm_br(dm_br.size());
  for(int i = 0; i < df_br.size(); i++){
    muf_br(i) = mx_mat_f(df_age_br(i)-1, df_time_br(i)-1) * exp(muf_br_tp(i)) * Ef_br(i);
  }
  for(int i = 0; i < dm_br.size(); i++){
    mum_br(i) = mx_mat_m(dm_age_br(i)-1, dm_time_br(i)-1) * exp(mum_br_tp(i)) * Em_br(i);
  }

  varf_br = muf_br * (1 + muf_br / exp(log_dispersion_br));
  varm_br = mum_br * (1 + mum_br / exp(log_dispersion_br));
  nll -= dnbinom2(df_br, muf_br, varf_br, 1).sum();
  nll -= dnbinom2(dm_br, mum_br, varm_br, 1).sum();

  matrix<Type> weights_mat_f(mx_mat_f.rows() - open_idx + 1, n_periods);
  weights_mat_f = Type(1.0) / (Type(1.0) + Type(0.5) * mx_mat_f.block(open_idx - 1, 0, mx_mat_f.rows() - open_idx + 1, n_periods).array());
  
  for(int i=0; i < n_periods; i++) {
  	for(int j=open_idx - 1; j < mx_mat_f.rows()-1; j++) {
		for(int k=j-open_idx+2; k < mx_mat_f.rows()-open_idx+1; k++) {
  			weights_mat_f(k,i) *= (Type(1.0) - Type(0.5) * mx_mat_f(j,i)) / (Type(1.0) + Type(0.5) * mx_mat_f(j,i));
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
  weights_mat_m = Type(1.0) / (Type(1.0) + Type(0.5) * mx_mat_m.block(open_idx - 1, 0, mx_mat_m.rows() - open_idx + 1, n_periods).array());
  
  for(int i=0; i < n_periods; i++) {
  	for(int j=open_idx - 1; j < mx_mat_m.rows()-1; j++) {
		for(int k=j-open_idx+2; k < mx_mat_m.rows()-open_idx+1; k++) {
  			weights_mat_m(k,i) *= (Type(1.0) - Type(0.5) * mx_mat_m(j,i)) / (Type(1.0) + Type(0.5) * mx_mat_m(j,i));
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
    sx_mat_f(0,i) = Type(1.0) / (Type(1.0) + Type(0.65) * mx_mat_f(0,i));
    sx_mat_m(0,i) = Type(1.0) / (Type(1.0) + Type(0.65) * mx_mat_m(0,i));

    sx_mat_f(1,i) = (Type(1.0) - Type(0.35) * mx_mat_f(0,i)) / (Type(1.0) + Type(0.5) * mx_mat_f(1,i));
    sx_mat_m(1,i) = (Type(1.0) - Type(0.35) * mx_mat_m(0,i)) / (Type(1.0) + Type(0.5) * mx_mat_m(1,i));
    //UDD open age group sx
    sx_mat_f(mx_mat_aggr_f.rows(), i) = (1 - Type(0.5) * mx_mat_aggr_f(mx_mat_aggr_f.rows() - 1, i)) / ((1 + Type(0.5) * mx_mat_aggr_f(mx_mat_aggr_f.rows() - 1, i)));
    sx_mat_m(mx_mat_aggr_m.rows(), i) = (1 - Type(0.5) * mx_mat_aggr_m(mx_mat_aggr_m.rows() - 1, i)) / ((1 + Type(0.5) * mx_mat_aggr_m(mx_mat_aggr_m.rows() - 1, i)));
    
    for(int j=1; j < (mx_mat_aggr_f.rows() - 1); j++) {
      sx_mat_f(j+1,i) = (Type(1.0) - Type(0.5) * mx_mat_f(j,i)) / (Type(1.0) + Type(0.5) * mx_mat_f(j+1,i));
      sx_mat_m(j+1,i) = (Type(1.0) - Type(0.5) * mx_mat_m(j,i)) / (Type(1.0) + Type(0.5) * mx_mat_m(j+1,i));
     }
   }

  Type sigma_logpop_m_base(exp(-0.5 * log_tau2_logpop(0)));
  Type sigma_logpop_m(exp(-0.5 * log_tau2_logpop(1)));
  Type sigma_logpop_f_base(exp(-0.5 * log_tau2_logpop(0)));
  Type sigma_logpop_f(exp(-0.5 * log_tau2_logpop(1)));

  //inverse gamma prior for variance with shape=1 and scale=0.0109  
  nll -= dlgamma(log_tau2_logpop, Type(1.0), Type(1.0 / 0.1), true).sum();

  nll -= dnorm(log_basepop_f, log_basepop_mean_f, sigma_logpop_f_base, true).sum();
  vector<Type> basepop_f(exp(log_basepop_f));
  nll -= dnorm(log_basepop_m, log_basepop_mean_m, sigma_logpop_m_base, true).sum();
  vector<Type> basepop_m(exp(log_basepop_m));

  int nnn = D_firstrow.size();
  // int nnn_ARIMA = D_firstrow_ARIMA.size();

  // prior for log(fx)
  Type rho_fx_age = invlogit(logit_rho_fx_age);
  Type rho_fx_time = invlogit(logit_rho_fx_time);

  nll -= AR2_rho_prior_base1(rho_fx_age, theta_rho_fx_age_base1);
  nll -= AR2_rho_prior_base1(rho_fx_time, theta_rho_fx_time_base1);

  nll -= Type(2.0) * log(rho_fx_age) - logit_rho_fx_age;
  nll -= Type(2.0) * log(rho_fx_time) - logit_rho_fx_time;
  
  Type theta_marginal_fx = -log(0.01) * sqrt(get_first_var_2d(D_firstrow , ARk(vector<Type>(rho_fx_time * rho_vec)).cov(nnn), ARk(vector<Type>(rho_fx_age * rho_vec)).cov(nnn))) / upper_marginal_sd_fx;
  Type lambda_fx = exp(log_marginal_lambda_fx);
  nll -= gumbel_density(lambda_fx, theta_marginal_fx);
  nll -= log_marginal_lambda_fx;
  Type sigma_fx = exp(-0.5 * log_marginal_lambda_fx);

  MapMatrixXXt fx_spline_mat(log_fx_spline_params.data(), log_fx_spline_params.size() / no_basis_time_int, no_basis_time_int);
  array<Type> fx_spline_array(fx_spline_mat.rows(), fx_spline_mat.cols());
  fx_spline_array = fx_spline_mat.array();
  nll += SCALE(SEPARABLE( ARk(vector<Type>(rho_fx_time * rho_vec)), ARk(vector<Type>(rho_fx_age * rho_vec)) ), sigma_fx)(fx_spline_array);

  vector<Type> fx(exp(log_fx_mean + D_agetime_fert * log_fx_spline_params));
  MapMatrixXXt fx_mat(fx.data(), n_fx, n_periods);

  // prior for gx
  Type rho_gx_age = invlogit(logit_rho_gx_age);
  Type rho_gx_time = invlogit(logit_rho_gx_time);

  nll -= AR2_rho_prior_base1(rho_gx_age, theta_rho_gx_age_base1);
  nll -= AR2_rho_prior_base1(rho_gx_time, theta_rho_gx_time_base1);

  nll -= Type(2.0) * log(rho_gx_age) - logit_rho_gx_age;
  nll -= Type(2.0) * log(rho_gx_time) - logit_rho_gx_time;

  Type theta_marginal_gx = -log(0.01) * sqrt(get_first_var_2d(D_firstrow , ARk(vector<Type>(rho_gx_time * rho_vec)).cov(nnn), ARk(vector<Type>(rho_gx_age * rho_vec)).cov(nnn))) / upper_marginal_sd_gx;
  Type lambda_gx = exp(log_marginal_lambda_gx);
  nll -= gumbel_density(lambda_gx, theta_marginal_gx);
  nll -= log_marginal_lambda_gx;
  Type sigma_gx = exp(-0.5 * log_marginal_lambda_gx);

  vector<Type> gx_f = D_agetime * gx_f_spline_params;
  vector<Type> gx_m = D_agetime * gx_m_spline_params;
  int basis_age = gx_f_spline_params.size() / no_basis_time_int;
  MapMatrixXXt gx_f_spline_mat(gx_f_spline_params.data(), basis_age, no_basis_time_int);
  MapMatrixXXt gx_m_spline_mat(gx_m_spline_params.data(), basis_age, no_basis_time_int);
  array<Type> gx_f_spline_array(gx_f_spline_mat.rows(), gx_f_spline_mat.cols());
  array<Type> gx_m_spline_array(gx_m_spline_mat.rows(), gx_m_spline_mat.cols());
  gx_f_spline_array = gx_f_spline_mat.array();
  gx_m_spline_array = gx_m_spline_mat.array();

  nll += SCALE(SEPARABLE( ARk(vector<Type>(rho_gx_time * rho_vec)), ARk(vector<Type>(rho_gx_age * rho_vec)) ), sigma_gx)(gx_f_spline_array);
  nll += SCALE(SEPARABLE( ARk(vector<Type>(rho_gx_time * rho_vec)), ARk(vector<Type>(rho_gx_age * rho_vec)) ), sigma_gx)(gx_m_spline_array);

  MapMatrixXXt gx_mat_f(gx_f.data(), basepop_f.size(), n_periods);
  MapMatrixXXt gx_mat_m(gx_m.data(), basepop_m.size(), n_periods);

  // population projection
  PopulationProjection<Type> proj(ccmpp<Type>(basepop_f, sx_mat_f, fx_mat, gx_mat_f,
					      srb, Type(1.0), fx_idx-1));
  PopulationProjection_m<Type> proj_m(ccmpp_m<Type>(basepop_m, sx_mat_m, gx_mat_m,
					      srb, Type(1.0), proj.births));

  matrix<Type> census_proj_mat_f(basepop_f.size(),census_year_idx.size());
  matrix<Type> census_proj_mat_m(basepop_m.size(),census_year_idx.size());
 
  //matrix<Type> census_proj_mat_f_oag(oag,census_year_idx.size());
  //matrix<Type> census_proj_mat_m_oag(oag,census_year_idx.size());

  for(int i = 0; i < census_year_idx.size(); i++) {	
	census_proj_mat_f.col(i) = vector<Type>(proj.population.col(census_year_idx[i] - 1));
		//.pow(Type(1.0) - census_year_grow_idx[i] / interval) * vector<Type>(proj.population.col(census_year_idx[i])).pow(census_year_grow_idx[i] / interval);
	census_proj_mat_m.col(i) = vector<Type>(proj_m.population.col(census_year_idx[i] - 1));
		//.pow(Type(1.0) - census_year_grow_idx[i] / interval) * vector<Type>(proj_m.population.col(census_year_idx[i])).pow(census_year_grow_idx[i] / interval);

	//census_proj_mat_f_oag.block(0, i, oag-1, 1) = census_proj_mat_f.block(0, i, oag-1, 1);
	//census_proj_mat_m_oag.block(0, i, oag-1, 1) = census_proj_mat_m.block(0, i, oag-1, 1);
	//census_proj_mat_f_oag(oag-1, i)  = census_proj_mat_f.col(i).segment(oag-1, basepop_f.size() - oag + 1).sum();
	//census_proj_mat_m_oag(oag-1, i)  = census_proj_mat_m.col(i).segment(oag-1, basepop_m.size() - oag + 1).sum();
  }
 
  matrix<Type> census_proj_mat_f_5(basepop_f.size(),census_year_idx_5.size());
  matrix<Type> census_proj_mat_m_5(basepop_m.size(),census_year_idx_5.size());
 
  //matrix<Type> census_proj_mat_f_oag_5(oag, census_year_idx_5.size());
  //matrix<Type> census_proj_mat_m_oag_5(oag, census_year_idx_5.size());

  matrix<Type> census_proj_mat_f_oag_5_aggr(n_agegrp_5, census_year_idx_5.size());
  matrix<Type> census_proj_mat_m_oag_5_aggr(n_agegrp_5, census_year_idx_5.size());

  for(int i = 0; i < census_year_idx_5.size(); i++) {	
	census_proj_mat_f_5.col(i) = vector<Type>(proj.population.col(census_year_idx_5[i] - 1));
	census_proj_mat_m_5.col(i) = vector<Type>(proj_m.population.col(census_year_idx_5[i] - 1));
	
	//census_proj_mat_f_oag_5.block(0, i, oag-1, 1) = census_proj_mat_f_5.block(0, i, oag-1, 1);
	//census_proj_mat_m_oag_5.block(0, i, oag-1, 1) = census_proj_mat_m_5.block(0, i, oag-1, 1);
	//census_proj_mat_f_oag_5(oag-1, i)  = census_proj_mat_f_5.col(i).segment(oag-1, basepop_f.size() - oag + 1).sum();
	//census_proj_mat_m_oag_5(oag-1, i)  = census_proj_mat_m_5.col(i).segment(oag-1, basepop_m.size() - oag + 1).sum();

	census_proj_mat_f_oag_5_aggr.col(i) = sumByAG1D(vector<Type>(census_proj_mat_f_5.col(i)), census_year_group_idx_5, n_agegrp_5);
	census_proj_mat_m_oag_5_aggr.col(i) = sumByAG1D(vector<Type>(census_proj_mat_m_5.col(i)), census_year_group_idx_5, n_agegrp_5);
  }

  // likelihood for log census counts
  for(int i = 0; i < census_year_idx.size(); i++) {
    nll -= dnorm(vector<Type>(census_log_pop_f.block(pop_start(i) - 1, i, pop_end(i) - pop_start(i) + 1, 1)),
  		 log(vector<Type>(census_proj_mat_f.block(pop_start(i) - 1, i, pop_end(i) - pop_start(i) + 1, 1))),
  		 sigma_logpop_f, true).sum();
    nll -= dnorm(vector<Type>(census_log_pop_m.block(pop_start(i) - 1, i, pop_end(i) - pop_start(i) + 1, 1)),
  		 log(vector<Type>(census_proj_mat_m.block(pop_start(i) - 1, i, pop_end(i) - pop_start(i) + 1, 1))),
  		 sigma_logpop_m, true).sum();
  }

  for(int i = 0; i < census_year_idx_5.size(); i++) {
    nll -= dnorm(vector<Type>(census_log_pop_f_5.block(pop_start_5(i) - 1, i, pop_end_5(i) - pop_start_5(i) + 1, 1)),
  		 log(vector<Type>(census_proj_mat_f_oag_5_aggr.block(pop_start_5(i) - 1, i, pop_end_5(i) - pop_start_5(i) + 1, 1))),
  		 sigma_logpop_f, true).sum();
    nll -= dnorm(vector<Type>(census_log_pop_m_5.block(pop_start_5(i) - 1, i, pop_end_5(i) - pop_start_5(i) + 1, 1)),
  		 log(vector<Type>(census_proj_mat_m_oag_5_aggr.block(pop_start_5(i) - 1, i, pop_end_5(i) - pop_start_5(i) + 1, 1))),
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

  // census deaths
  nll -= dnorm(log_dispersion_census_deaths, Type(0.0), Type(5.0), 1).sum();
  nll -= dnorm(log_dispersion_census_deaths_5, Type(0.0), Type(5.0), 1).sum();

  matrix<Type> census_deaths_mat_f(basepop_f.size(),census_deaths_year_idx.size());
  matrix<Type> census_deaths_mat_m(basepop_m.size(),census_deaths_year_idx.size());
  matrix<Type> census_deaths_varmat_f(basepop_f.size(),census_deaths_year_idx.size());
  matrix<Type> census_deaths_varmat_m(basepop_m.size(),census_deaths_year_idx.size());

  for(int i = 0; i < census_deaths_year_idx.size(); i++) {
	census_deaths_mat_f.col(i) = vector<Type>(proj_period_deaths_f.col(census_deaths_year_idx[i] - 2));
	census_deaths_mat_m.col(i) = vector<Type>(proj_period_deaths_m.col(census_deaths_year_idx[i] - 2));

	census_deaths_varmat_f.col(i) = vector<Type>(census_deaths_mat_f.col(i)) * (Type(1.0) + vector<Type>(census_deaths_mat_f.col(i)) / exp(log_dispersion_census_deaths(1)));
	census_deaths_varmat_m.col(i) = vector<Type>(census_deaths_mat_m.col(i)) * (Type(1.0) + vector<Type>(census_deaths_mat_m.col(i)) / exp(log_dispersion_census_deaths(0)));

    nll -= dnbinom2(vector<Type>(census_deaths_f.block(deaths_start(i) - 1, i, deaths_end(i) - deaths_start(i) + 1, 1)),
  		 vector<Type>(census_deaths_mat_f.block(deaths_start(i) - 1, i, deaths_end(i) - deaths_start(i) + 1, 1)),
  		 vector<Type>(census_deaths_varmat_f.block(deaths_start(i) - 1, i, deaths_end(i) - deaths_start(i) + 1, 1)),
                 1).sum();
    nll -= dnbinom2(vector<Type>(census_deaths_m.block(deaths_start(i) - 1, i, deaths_end(i) - deaths_start(i) + 1, 1)),
  		 vector<Type>(census_deaths_mat_m.block(deaths_start(i) - 1, i, deaths_end(i) - deaths_start(i) + 1, 1)),
  		 vector<Type>(census_deaths_varmat_m.block(deaths_start(i) - 1, i, deaths_end(i) - deaths_start(i) + 1, 1)),
                 1).sum();

//    nll -= dpois(vector<Type>(census_deaths_f.block(deaths_start(i) - 1, i, deaths_end(i) - deaths_start(i) + 1, 1)), 
//		vector<Type>(census_deaths_mat_f.block(deaths_start(i) - 1, i, deaths_end(i) - deaths_start(i) + 1, 1)), 1).sum();
//    nll -= dpois(vector<Type>(census_deaths_m.block(deaths_start(i) - 1, i, deaths_end(i) - deaths_start(i) + 1, 1)), 
//		vector<Type>(census_deaths_mat_m.block(deaths_start(i) - 1, i, deaths_end(i) - deaths_start(i) + 1, 1)), 1).sum();
  }
 
  matrix<Type> census_deaths_mat_f_5(basepop_f.size(),census_deaths_year_idx_5.size());
  matrix<Type> census_deaths_mat_m_5(basepop_m.size(),census_deaths_year_idx_5.size()); 
  matrix<Type> census_deaths_mat_f_5_aggr(n_agegrp_5_deaths, census_deaths_year_idx_5.size());
  matrix<Type> census_deaths_mat_m_5_aggr(n_agegrp_5_deaths, census_deaths_year_idx_5.size());
  matrix<Type> census_deaths_varmat_f_5_aggr(n_agegrp_5_deaths, census_deaths_year_idx_5.size());
  matrix<Type> census_deaths_varmat_m_5_aggr(n_agegrp_5_deaths, census_deaths_year_idx_5.size());

  for(int i = 0; i < census_deaths_year_idx_5.size(); i++) {	
	census_deaths_mat_f_5.col(i) = vector<Type>(proj_period_deaths_f.col(census_deaths_year_idx_5[i] - 2));
	census_deaths_mat_m_5.col(i) = vector<Type>(proj_period_deaths_m.col(census_deaths_year_idx_5[i] - 2));
	census_deaths_mat_f_5_aggr.col(i) = sumByAG1D(vector<Type>(census_deaths_mat_f_5.col(i)), census_deaths_group_idx_5, n_agegrp_5_deaths);
	census_deaths_mat_m_5_aggr.col(i) = sumByAG1D(vector<Type>(census_deaths_mat_m_5.col(i)), census_deaths_group_idx_5, n_agegrp_5_deaths);

        census_deaths_varmat_f_5_aggr.col(i) = vector<Type>(census_deaths_mat_f_5_aggr.col(i)) * (Type(1.0) + vector<Type>(census_deaths_mat_f_5_aggr.col(i)) / exp(log_dispersion_census_deaths_5(1)));
	census_deaths_varmat_m_5_aggr.col(i) = vector<Type>(census_deaths_mat_m_5_aggr.col(i)) * (Type(1.0) + vector<Type>(census_deaths_mat_m_5_aggr.col(i)) / exp(log_dispersion_census_deaths_5(0)));

    nll -= dnbinom2(vector<Type>(census_deaths_f_5.block(deaths_start_5(i) - 1, i, deaths_end_5(i) - deaths_start_5(i) + 1, 1)),
  		 vector<Type>(census_deaths_mat_f_5_aggr.block(deaths_start_5(i) - 1, i, deaths_end_5(i) - deaths_start_5(i) + 1, 1)),
  		 vector<Type>(census_deaths_varmat_f_5_aggr.block(deaths_start_5(i) - 1, i, deaths_end_5(i) - deaths_start_5(i) + 1, 1)),
                 1).sum();
    nll -= dnbinom2(vector<Type>(census_deaths_m_5.block(deaths_start_5(i) - 1, i, deaths_end_5(i) - deaths_start_5(i) + 1, 1)),
  		 vector<Type>(census_deaths_mat_m_5_aggr.block(deaths_start_5(i) - 1, i, deaths_end_5(i) - deaths_start_5(i) + 1, 1)),
  		 vector<Type>(census_deaths_varmat_m_5_aggr.block(deaths_start_5(i) - 1, i, deaths_end_5(i) - deaths_start_5(i) + 1, 1)),
                 1).sum();
  }

   
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

    REPORT(log_par_f);
	REPORT(log_par_m);
	
    REPORT(phi_f);
    REPORT(disp_f);
    REPORT(psi_f);
    REPORT(lambda_f);
    REPORT(delta_f);
    REPORT(epsilon_f);
    REPORT(A_f);
    REPORT(B_f);

    REPORT(phi_m);
    REPORT(disp_m);
    REPORT(psi_m);
    REPORT(lambda_m);
    REPORT(delta_m);
    REPORT(epsilon_m);
    REPORT(A_m);
    REPORT(B_m);

    REPORT(tp_params);
    REPORT(tp_params_br_age0);
    REPORT(tp_params_br_age1);
    REPORT(tp_params_br_age5);
    REPORT(tp_params_br_age10);

  return Type(nll);
}

