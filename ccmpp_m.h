#ifndef CCMPP_M_H
#define CCMPP_M_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

template <typename Type>
class PopulationProjection_m {

  typedef Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> MatrixXXT;
  typedef Eigen::Array<Type, Eigen::Dynamic, Eigen::Dynamic> ArrayXXT;
  typedef Eigen::Array<Type, 1, Eigen::Dynamic> Array1XT;

 public: 

  const int n_ages;
  const int n_periods;
  const Type interval;

  const ArrayXXT sx;
  const ArrayXXT gx;
  const Array1XT srb;
  const MatrixXXT births;  

  MatrixXXT population;
  MatrixXXT cohort_deaths;
  Eigen::Matrix<Type, 1, Eigen::Dynamic> infants;
  MatrixXXT migrations;
    
 PopulationProjection_m(const int n_ages,
		      const int n_periods,
		      const Type interval,
		      const Eigen::Matrix<Type, Eigen::Dynamic, 1>& basepop,
		      const Eigen::Array<Type, Eigen::Dynamic, Eigen::Dynamic>& sx,
		      const Eigen::Array<Type, Eigen::Dynamic, Eigen::Dynamic>& gx,
		      const Eigen::Array<Type, Eigen::Dynamic, 1>& srb,
		      const Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>& births) :
    n_ages{ n_ages },
    n_periods{ n_periods },
    interval{ interval },
    sx{ sx },
    gx{ gx },
    srb{ srb },
    births{ births },
    population{ MatrixXXT(n_ages, n_periods + 1) },
    cohort_deaths{ MatrixXXT(n_ages+1, n_periods) },
    infants{ Eigen::Matrix<Type, 1, Eigen::Dynamic>(n_periods) },
    migrations{ MatrixXXT(n_ages, n_periods) }
    {
      population.col(0) = basepop;
    };

    void step_projection(int step);

    MatrixXXT period_deaths();
  
};

template <typename Type>
void PopulationProjection_m<Type>::step_projection(int step) {

  typedef Eigen::Map<Eigen::Array<Type, Eigen::Dynamic, 1> > MapArrayXT;
  typedef Eigen::Map<const Eigen::Array<Type, Eigen::Dynamic, 1> > MapConstArrayXT;

  MapConstArrayXT sx_t(sx.col(step).data(), n_ages + 1);
  MapConstArrayXT gx_t(gx.col(step).data(), n_ages);
  MapConstArrayXT births_t(births.col(step).data(),births.rows());

  MapArrayXT population_t(population.col(step+1).data(), population.rows());
  MapArrayXT migrations_t(migrations.col(step).data(), migrations.rows());
  MapArrayXT cohort_deaths_t(cohort_deaths.col(step).data(), cohort_deaths.rows());
  MapArrayXT infants_t(infants.col(step).data(), infants.rows());

  population_t = population.col(step);
  
  migrations_t = population_t * gx_t;
  population_t += 0.5 * migrations_t;
  
  cohort_deaths_t.segment(1, n_ages) = population_t * (1.0 - sx_t.segment(1, n_ages));

  Type open_age_survivors = population_t(n_ages-1) - cohort_deaths_t(n_ages);
  for(int age = n_ages-1; age > 0; age--) {
    population_t(age) = population_t(age-1) - cohort_deaths_t(age);
  }
  population_t(n_ages-1) += open_age_survivors;
  
  infants_t = births_t.sum() * srb(step) / (1.0 + srb(step));
  cohort_deaths_t(0) = infants_t(0) * (1.0 - sx_t(0));
  population_t(0) = infants_t(0) - cohort_deaths_t(0);
  
  population_t += 0.5 * migrations_t;
}

template <typename Type>
Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> PopulationProjection_m<Type>::period_deaths() {

  typedef Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> MatrixXXT;
  
  MatrixXXT period_deaths(0.5 * cohort_deaths.topRows(n_ages));
  period_deaths += 0.5 * cohort_deaths.bottomRows(n_ages);
  period_deaths.row(0) += 0.5 * cohort_deaths.row(0);
  period_deaths.row(n_ages-1) += 0.5 * cohort_deaths.row(n_ages);

  return period_deaths;
}

  
template <typename Type>
PopulationProjection_m<Type>
ccmpp_m(const Eigen::Matrix<Type, Eigen::Dynamic, 1>& basepop,
      const Eigen::Array<Type, Eigen::Dynamic, Eigen::Dynamic>& sx,
      const Eigen::Array<Type, Eigen::Dynamic, Eigen::Dynamic>& gx,
      const Eigen::Array<Type, Eigen::Dynamic, 1>& srb,
      const Type interval,
      const Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>& births) {

  int n_periods(sx.cols());
  int n_ages(basepop.rows());

  PopulationProjection_m<Type> proj(n_ages, n_periods, interval,
				  basepop, sx, gx, srb, births);

  for(int step = 0; step < n_periods; step++) {
    proj.step_projection(step);
  }

  return proj;
}

#endif
