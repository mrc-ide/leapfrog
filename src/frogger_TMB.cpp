#include <Rcpp.h>
#include <TMB.hpp>

#include "frogger.hpp"
#include "intermediate_data.hpp"
#include "model_variants.hpp"
#include "state_space.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj
template<typename ModelVariant, typename Type>
int simulate_model_TMB(objective_function<Type>* obj) {

  DATA_INTEGER(proj_years)
  DATA_INTEGER(hiv_steps)
  DATA_IVECTOR(save_steps)

  DATA_INTEGER(t_ART_start)

  constexpr auto ss = leapfrog::StateSpace<ModelVariant>();

  const leapfrog::Options<Type> opts = {
      hiv_steps,
      t_ART_start - 1,
      ss.base.hAG
  };

  // setup params

  // run model and return state, for now returning 0
  return 0;
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_STRING(model_variant)
  
  PARAMETER(x)

  if (model_variant == "ChildModel") {
    auto output_state = simulate_model_TMB<leapfrog::ChildModel, Type>(this);
  } else if (model_variant == "BaseModelFullAgeStratification") {
    auto output_state = simulate_model_TMB<leapfrog::BaseModelFullAgeStratification, Type>(this);
  } else if (model_variant == "BaseModelCoarseAgeStratification") {
    auto output_state = simulate_model_TMB<leapfrog::BaseModelCoarseAgeStratification, Type>(this);
  } else {
    throw std::invalid_argument("Invalid model variant " + model_variant + " must be one of " +
                                "'BaseModelFullAgeStratification', 'BaseModelCoarseAgeStratification' or 'ChildModel'");
  }

  return 0;
}
