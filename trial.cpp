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

matrix<Type> A;
matrix<Type> B;
matrix<Type> C;
vector<Type> d;
vector<Type> e;

MapMatrixXXt M(d.data(), 2, 2);

MVNORM_t<Type> AAA(M);

Type nll = AAA(A.col(0));
}

