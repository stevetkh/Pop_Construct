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


template<class Type>
Type objective_function<Type>::operator() ()
{
 
DATA_VECTOR(aa);
DATA_IVECTOR(bb);
DATA_INTEGER(c);

Type nll=0.0;

vector<Type> haha = sumByAG1D(aa, bb, c);

return Type(nll);
}


