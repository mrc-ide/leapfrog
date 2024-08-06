#include <Rcpp.h>
#include <TMB.hpp>

#include "frogger.hpp"
#include "intermediate_data.hpp"
#include "model_variants.hpp"
#include "state_space.hpp"
#include "generated/model_input_TMB_new.hpp"

template<typename ModelVariant>
struct TMBDataStruct {
  template<typename Type>
  using type = TMBData<Type, ModelVariant>;
};

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
template<typename ModelVariant, typename Type>
Type simulate_model_TMB(objective_function<Type>* obj) {

  DATA_INTEGER(proj_years)
  DATA_INTEGER(hiv_steps)
  DATA_IVECTOR(save_steps)
  DATA_INTEGER(t_ART_start)

  DATA_STRUCT(data_vars, typename TMBDataStruct<ModelVariant>::type);

  constexpr leapfrog::StateSpace<ModelVariant> ss = leapfrog::StateSpace<ModelVariant>();
  const leapfrog::Options<Type> opts = {
      hiv_steps,
      t_ART_start - 1,
      ss.base.hAG
  };

  leapfrog::Parameters<ModelVariant, Type> params = setup_tmb_params(data_vars, opts, proj_years);

  // PARAMETER_ARRAY(incidinput);
  PARAMETER(incidinput_scalar);
  Eigen::Tensor<Type, 1> incidinput_const(proj_years);
  Type input_clean = 0;
  Type bounded_incidinput = invlogit(incidinput_scalar);
  incidinput_const.setConstant(bounded_incidinput);
  params.base.incidence.total_rate = Eigen::TensorMap<const Eigen::Tensor<Type, 1>>(incidinput_const.data(), proj_years);
  DATA_ARRAY(basepop);
  Eigen::TensorMap<const Eigen::Tensor<Type, 3>> base_pop = Eigen::TensorMap<const Eigen::Tensor<Type, 3>>(basepop.data(), ss.base.pAG, ss.base.NS, proj_years);

  DATA_ARRAY(incidinput_data)

  leapfrog::OutputState<ModelVariant, Type> state = leapfrog::run_model<ModelVariant, Type>(proj_years, save_steps, params);

  Type nll = 0;
  DATA_SCALAR(sd);

  // for (int y = 0; y < proj_years; y++) {
    Type total_pop = 0;
    Type total_p_new_infections = 0;
    for (int a = 0; a < ss.base.pAG; a++) {
      for (int g = 0; g < ss.base.NS; g++) {
        total_pop += base_pop(a, g, 30);
        total_p_new_infections += state.base.p_infections(a, g, 30);
      }
    }
    Type total_new_infections = total_pop * incidinput_data(30);
    Rcout << "Proj year: " << 30 << std::endl;
    Rcout << "New infections: " << total_new_infections << std::endl;
    Rcout << "Projected new infections: " << total_p_new_infections << std::endl;
    Rcout << "I am testing: " << params.base.incidence.total_rate(30) << std::endl;
    Rcout << std::endl;

    nll += -dnorm(total_new_infections, total_p_new_infections, sd, true);
  // }

  return nll;
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_STRING(model_variant)
  
  if (model_variant == "ChildModel") {
    return simulate_model_TMB<leapfrog::ChildModel, Type>(this);
  } else if (model_variant == "BaseModelFullAgeStratification") {
    return simulate_model_TMB<leapfrog::BaseModelFullAgeStratification, Type>(this);
  } else if (model_variant == "BaseModelCoarseAgeStratification") {
    return simulate_model_TMB<leapfrog::BaseModelCoarseAgeStratification, Type>(this);
  } else {
    throw std::invalid_argument("Invalid model variant " + model_variant + " must be one of " +
                                "'BaseModelFullAgeStratification', 'BaseModelCoarseAgeStratification' or 'ChildModel'");
  }

  return 0;
}
