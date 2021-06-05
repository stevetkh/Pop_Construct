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

  PARAMETER_VECTOR(log_tau2_logpop); //1 = WPP males, 2 = data males, 3 = WPP females, 4 = data females
  PARAMETER(log_lambda_fx);
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
  DATA_VECTOR(thiele_age);

  DATA_VECTOR(log_phi_mean_f);
  DATA_VECTOR(log_psi_mean_f);
  DATA_VECTOR(log_A_mean_f);
  DATA_VECTOR(log_B_mean_f);
  DATA_SCALAR(log_lambda_mean_f);
  DATA_SCALAR(log_delta_mean_f);
  DATA_SCALAR(log_epsilon_mean_f);

  DATA_VECTOR(log_phi_mean_m);
  DATA_VECTOR(log_psi_mean_m);
  DATA_VECTOR(log_A_mean_m);
  DATA_VECTOR(log_B_mean_m);
  DATA_SCALAR(log_lambda_mean_m);
  DATA_SCALAR(log_delta_mean_m);
  DATA_SCALAR(log_epsilon_mean_m);

  DATA_SPARSE_MATRIX(full_penal_time);

  DATA_SPARSE_MATRIX(full_penal_gx);

  DATA_SPARSE_MATRIX(full_penal_fx);

  DATA_SPARSE_MATRIX(D_time);
  DATA_SPARSE_MATRIX(D_agetime);
  DATA_SPARSE_MATRIX(D_agetime_fert);

  DATA_SPARSE_MATRIX(penal_tp);
  DATA_SPARSE_MATRIX(penal_tp_0);
  DATA_SPARSE_MATRIX(null_penal_tp);
  DATA_VECTOR(tp_mean);
  
  PARAMETER_VECTOR(log_dispersion); //1 = male, 2 = female

  PARAMETER(log_lambda_tp);
  PARAMETER_VECTOR(tp_params);
  PARAMETER(tp_slope);
  PARAMETER(tp_params_5);
  PARAMETER(tp_params_10);

  PARAMETER_VECTOR(log_phi_f_spline_params); 
  PARAMETER_VECTOR(log_psi_f_spline_params);
  PARAMETER_VECTOR(log_lambda_f_spline_params);
  PARAMETER_VECTOR(log_delta_f_spline_params);
  PARAMETER_VECTOR(log_epsilon_f_spline_params);
  PARAMETER_VECTOR(log_A_f_spline_params);
  PARAMETER_VECTOR(log_B_f_spline_params);

  PARAMETER_VECTOR(log_phi_m_spline_params); 
  PARAMETER_VECTOR(log_psi_m_spline_params);
  PARAMETER_VECTOR(log_lambda_m_spline_params);
  PARAMETER_VECTOR(log_delta_m_spline_params);
  PARAMETER_VECTOR(log_epsilon_m_spline_params);
  PARAMETER_VECTOR(log_A_m_spline_params);
  PARAMETER_VECTOR(log_B_m_spline_params);

  PARAMETER(log_lambda_phi);
  PARAMETER(log_lambda_psi);
  PARAMETER(log_lambda_A);
  PARAMETER(log_lambda_B);

  PARAMETER(log_lambda_lambda);
  PARAMETER(log_lambda_delta);
  PARAMETER(log_lambda_epsilon);

  DATA_SCALAR(theta_fx);
  DATA_SCALAR(theta_marginal_gx);
  DATA_SCALAR(theta_phi);
  DATA_SCALAR(theta_psi);
  DATA_SCALAR(theta_lambda);
  DATA_SCALAR(theta_delta);
  DATA_SCALAR(theta_epsilon);
  DATA_SCALAR(theta_A);
  DATA_SCALAR(theta_B);
  
  DATA_SCALAR(theta_tp);

  Type nll(0.0);

  //inverse gamma prior for variance with shape=1 and scale=0.0109  
  nll -= dlgamma(log_tau2_logpop, Type(1.0), Type(1.0 / 0.1), true).sum();
  Type sigma_logpop_m_base(exp(-0.5 * log_tau2_logpop(0)));
  Type sigma_logpop_m(exp(-0.5 * log_tau2_logpop(1)));
  Type sigma_logpop_f_base(exp(-0.5 * log_tau2_logpop(2)));
  Type sigma_logpop_f(exp(-0.5 * log_tau2_logpop(3)));

  Type lambda_fx = exp(log_lambda_fx);
  nll -= gumbel_density(lambda_fx, theta_fx);
  nll -= log_lambda_fx;

  Type marginal_lambda_gx = exp(log_marginal_lambda_gx);
  nll -= gumbel_density(marginal_lambda_gx, theta_marginal_gx);
  nll -= log_marginal_lambda_gx;
  Type sigma_gx = exp(-0.5 * log_marginal_lambda_gx);

  nll -= dnorm(log_basepop_f, log_basepop_mean_f, sigma_logpop_f_base, true).sum();
  vector<Type> basepop_f(exp(log_basepop_f));
  nll -= dnorm(log_basepop_m, log_basepop_mean_m, sigma_logpop_m_base, true).sum();
  vector<Type> basepop_m(exp(log_basepop_m));

  Type lambda_phi = exp(log_lambda_phi);
  Type lambda_psi = exp(log_lambda_psi);
  Type lambda_A = exp(log_lambda_A);
  Type lambda_B = exp(log_lambda_B);

  nll -= gumbel_density(lambda_phi, theta_phi);
  nll -= gumbel_density(lambda_psi, theta_psi);
  nll -= gumbel_density(lambda_A, theta_A);
  nll -= gumbel_density(lambda_B, theta_B);
  nll -= log_lambda_phi;
  nll -= log_lambda_psi;
  nll -= log_lambda_A;
  nll -= log_lambda_B;

  Type lambda_lambda = exp(log_lambda_lambda);
  Type lambda_delta = exp(log_lambda_delta);
  Type lambda_epsilon = exp(log_lambda_epsilon);

  nll -= gumbel_density(lambda_lambda, theta_lambda);
  nll -= gumbel_density(lambda_delta, theta_delta);
  nll -= gumbel_density(lambda_epsilon, theta_epsilon);
  nll -= log_lambda_lambda;
  nll -= log_lambda_delta;
  nll -= log_lambda_epsilon;

  Type lambda_tp = exp(log_lambda_tp);
  nll -= gumbel_density(lambda_tp, theta_tp);

  SparseMatrix<Type> QQ_tp = lambda_tp * penal_tp + 4 * penal_tp_0 + null_penal_tp;
  nll += GMRF(QQ_tp)(tp_params - tp_slope * tp_mean);
  nll -= dnorm(tp_slope, Type(0.0), Type(0.3), true);
  nll -= dnorm(tp_params_5, Type(0.0), Type(0.3), true);
  nll -= dnorm(tp_params_10, Type(0.0), Type(0.3), true);
  tp_params(5) += tp_params_5;
  tp_params(10) += tp_params_10;

  nll -= dnorm(log_dispersion, Type(0.0), Type(5.0), 1).sum();

  SparseMatrix<Type> QQ_phi =  lambda_phi * full_penal_time;
  nll += GMRF(QQ_phi)(log_phi_f_spline_params) + GMRF(QQ_phi)(log_phi_m_spline_params);
  SparseMatrix<Type> QQ_psi =  lambda_psi * full_penal_time;
  nll += GMRF(QQ_psi)(log_psi_f_spline_params) + GMRF(QQ_psi)(log_psi_m_spline_params);
  SparseMatrix<Type> QQ_lambda =  lambda_lambda * full_penal_time;
  nll += GMRF(QQ_lambda)(log_lambda_f_spline_params) + GMRF(QQ_lambda)(log_lambda_m_spline_params);
  SparseMatrix<Type> QQ_delta =  lambda_delta * full_penal_time;
  nll += GMRF(QQ_delta)(log_delta_f_spline_params) + GMRF(QQ_delta)(log_delta_m_spline_params);
  SparseMatrix<Type> QQ_epsilon =  lambda_epsilon * full_penal_time;
  nll += GMRF(QQ_epsilon)(log_epsilon_f_spline_params) + GMRF(QQ_epsilon)(log_epsilon_m_spline_params);
  SparseMatrix<Type> QQ_A =  lambda_A * full_penal_time;
  nll += GMRF(QQ_A)(log_A_f_spline_params) + GMRF(QQ_A)(log_A_m_spline_params);
  SparseMatrix<Type> QQ_B =  lambda_B * full_penal_time;
  nll += GMRF(QQ_B)(log_B_f_spline_params) + GMRF(QQ_B)(log_B_m_spline_params);

  vector<Type> phi_f = exp(log_phi_mean_f + D_time * log_phi_f_spline_params);
  vector<Type> psi_f = exp(log_psi_mean_f + D_time * log_psi_f_spline_params);
  vector<Type> lambda_f = exp(log_lambda_mean_f + D_time * log_lambda_f_spline_params);
  vector<Type> delta_f = exp(log_delta_mean_f + D_time * log_delta_f_spline_params);
  vector<Type> epsilon_f = exp(log_epsilon_mean_f + D_time * log_epsilon_f_spline_params);
  vector<Type> A_f = exp(log_A_mean_f + D_time * log_A_f_spline_params);
  vector<Type> B_f = exp(log_B_mean_f + D_time * log_B_f_spline_params);

  vector<Type> phi_m = exp(log_phi_mean_m + D_time * log_phi_m_spline_params);
  vector<Type> psi_m = exp(log_psi_mean_m + D_time * log_psi_m_spline_params);
  vector<Type> lambda_m = exp(log_lambda_mean_m + D_time * log_lambda_m_spline_params);
  vector<Type> delta_m = exp(log_delta_mean_m + D_time * log_delta_m_spline_params);
  vector<Type> epsilon_m = exp(log_epsilon_mean_m + D_time * log_epsilon_m_spline_params);
  vector<Type> A_m = exp(log_A_mean_m + D_time * log_A_m_spline_params);
  vector<Type> B_m = exp(log_B_mean_m + D_time * log_B_m_spline_params);

  matrix<Type> mx_mat_f(thiele_age.size(), phi_f.size()); 
  matrix<Type> mx_mat_m(thiele_age.size(), phi_f.size()); 
  for(int i = 0; i < phi_f.size(); i++){
    mx_mat_f.col(i) = phi_f(i)*exp(-psi_f(i)*(thiele_age - Type(2.0))) + lambda_f(i)*exp(-delta_f(i)*((log(thiele_age)-log(epsilon_f(i)))*(log(thiele_age)-log(epsilon_f(i))))) + A_f(i)*exp(B_f(i)*(thiele_age - Type(92.0)));
    mx_mat_m.col(i) = phi_m(i)*exp(-psi_m(i)*(thiele_age - Type(2.0))) + lambda_m(i)*exp(-delta_m(i)*((log(thiele_age)-log(epsilon_m(i)))*(log(thiele_age)-log(epsilon_m(i))))) + A_m(i)*exp(B_m(i)*(thiele_age - Type(92.0)));
  }

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

  varf = muf * (1 + muf / exp(log_dispersion(1)));
  varm = mum * (1 + mum / exp(log_dispersion(0)));
  nll -= dnbinom2(df, muf, varf, 1).sum();
  nll -= dnbinom2(dm, mum, varm, 1).sum();

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
    sx_mat_f(0,i) = (Type(1.0) - Type(0.15) * mx_mat_f(0,i)) / (Type(1.0) + Type(0.5) * mx_mat_f(0,i));
    sx_mat_m(0,i) = (Type(1.0) - Type(0.15) * mx_mat_m(0,i)) / (Type(1.0) + Type(0.5) * mx_mat_m(0,i));
    //UDD open age group sx
    sx_mat_f(mx_mat_aggr_f.rows(), i) = (1 - Type(0.5) * mx_mat_aggr_f(mx_mat_aggr_f.rows() - 1, i)) / ((1 + Type(0.5) * mx_mat_aggr_f(mx_mat_aggr_f.rows() - 1, i)));
    sx_mat_m(mx_mat_aggr_m.rows(), i) = (1 - Type(0.5) * mx_mat_aggr_m(mx_mat_aggr_m.rows() - 1, i)) / ((1 + Type(0.5) * mx_mat_aggr_m(mx_mat_aggr_m.rows() - 1, i)));
    
    for(int j=0; j < (mx_mat_aggr_f.rows() - 1); j++) {
      sx_mat_f(j+1,i) = (Type(1.0) - Type(0.5) * mx_mat_f(j,i)) / (Type(1.0) + Type(0.5) * mx_mat_f(j+1,i));
      sx_mat_m(j+1,i) = (Type(1.0) - Type(0.5) * mx_mat_m(j,i)) / (Type(1.0) + Type(0.5) * mx_mat_m(j+1,i));
     }
   }

  // prior for log(fx)
  SparseMatrix<Type> QQ_fx = lambda_fx * full_penal_fx;
  nll += GMRF(QQ_fx)(log_fx_spline_params);
  vector<Type> fx(exp(log_fx_mean + D_agetime_fert * log_fx_spline_params));
  MapMatrixXXt fx_mat(fx.data(), n_fx, n_periods);

  // prior for gx
  vector<Type> gx_f = D_agetime * gx_f_spline_params;
  vector<Type> gx_m = D_agetime * gx_m_spline_params;
  int basis_age = gx_f_spline_params.size() / log_phi_f_spline_params.size();
  MapMatrixXXt gx_f_spline_mat(gx_f_spline_params.data(), basis_age, log_phi_f_spline_params.size());
  MapMatrixXXt gx_m_spline_mat(gx_m_spline_params.data(), basis_age, log_phi_m_spline_params.size());
  array<Type> gx_f_spline_array(gx_f_spline_mat.rows(), gx_f_spline_mat.cols());
  array<Type> gx_m_spline_array(gx_m_spline_mat.rows(), gx_m_spline_mat.cols());
  gx_f_spline_array = gx_f_spline_mat.array();
  gx_m_spline_array = gx_m_spline_mat.array();
  vector<Type> rho_vec(2);
  rho_vec(0) = 2 * 0.9;
  rho_vec(1) = -0.9;

  nll += SCALE(SEPARABLE( ARk(rho_vec), ARk(rho_vec) ), sigma_gx)(gx_f_spline_array);
  nll += SCALE(SEPARABLE( ARk(rho_vec), ARk(rho_vec) ), sigma_gx)(gx_m_spline_array);

  //SparseMatrix<Type> QQ_gx = lambda_gx * full_penal_gx;
  //nll += GMRF(QQ_gx)(gx_f_spline_params) + GMRF(QQ_gx)(gx_m_spline_params);

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
