#include <TMB.hpp>                                

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
Type d_gamma_innov(const Type& gamma, const Type& n)
{
  Type ll = 0.0;
  
  ll += (gamma - Type(1.0)) * (n - Type(1.5)) - (n - Type(3.0)) * log(gamma);
  ll += - log(Type(0.25)*gamma + Type(0.75)) - log(Type(0.5)*gamma + Type(0.5)) - log(Type(0.75)*gamma + Type(0.25));
  
  return sqrt(ll);
}


template <typename Type>
Type logdd_gamma_innov(const Type& gamma, const Type& n)
{
  Type ll = 0.0; 
  Type chunk = 0.0;
  Type chunk_abs = 0.0;

  chunk += n - Type(1.5) - (n - Type(3.0)) / gamma;
  chunk += - Type(0.25)/(Type(0.25)*gamma + Type(0.75)) - Type(0.5)/(Type(0.5)*gamma + Type(0.5)) - Type(0.75)/(Type(0.75)*gamma + Type(0.25));

  chunk_abs = sqrt(chunk * chunk); 


  ll += log(Type(0.5)) - log(d_gamma_innov(gamma, n)) + log(chunk_abs);

  return ll;
}


template <typename Type>
Type AR2_innov_prec(const Type& rho, int n)
{
  matrix<Type> I(n, n);
  Type tau2 = 0.0;
  Type det_init = 0.0;
  I.setZero();

  tau2 = (Type(1.0) + rho) / ((Type(1.0) - rho) * (Type(1.0) - rho) * (Type(1.0) + Type(3.0) * rho));

  det_init = Type(1.0) - 4 * rho * rho / ((Type(1.0) + rho) * (Type(1.0) + rho));
  
  I(0,0) = Type(1.0) / det_init;
  I(1,1) = Type(1.0) / det_init;
  I(0,1) = -Type(2.0) * rho / (det_init * (Type(1.0) + rho));
  I(1,0) = -Type(2.0) * rho / (det_init * (Type(1.0) + rho));

  for(int i = 2; i < n; i++){
    I(i,i) = tau2;
  }
  return I;
}

template <typename Type>
Type AR2_T(const Type& rho, int n)
{
  matrix<Type> I(n, n);
  I.setZero();

  I(0,0) = Type(1.0);
  I(1,1) = Type(1.0);  

  for(int i = 2; i < n; i++){
    I(i, i) = Type(1.0);
    I(i, i-1) = -Type(2.) * rho;
    I(i, i-2) = rho;
  }
   return I;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
using namespace Eigen;
using namespace density;

DATA_MATRIX(aa);
DATA_MATRIX(bb);

Type nll=0.0;

SparseMatrix<Type> F = asSparseMatrix( bb.transpose() * aa * bb);

return Type(nll);
}


