#pragma once

#include "types.hpp"

namespace leapfrog {

namespace internal {


template<typename ModelVariant, typename real_type>
void hc_adjust_art_initiates_for_mort_age_spec(int time_step,
                         const Parameters<ModelVariant, real_type> &pars,
                         const State<ModelVariant, real_type> &state_curr,
                         State<ModelVariant, real_type> &state_next,
                         IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "run_hiv_child_infections can only be called for model variants where run_child_model is true");
  constexpr auto ss = StateSpace<ModelVariant>().base;
  constexpr auto hc_ss = StateSpace<ModelVariant>().children;
  const auto cpars = pars.children.children;

  for (int s = 0; s <ss.NS; ++s) {
    for (int a = 0; a < pars.base.options.p_idx_fertility_first; ++a) {
      for (int hd = 0; hd < hc_ss.hc1DS; ++hd) {
        for (int cat = 0; cat < hc_ss.hcTT; ++cat) {
          state_next.children.hc_art_num += intermediate.children.hc_art_need(hd, cat, a, s);
        } // end hcTT
      } // end hc_ss.hc1DS
    } // end a
  } // end ss.NS

  //how many should initialize ART
  for (int s = 0; s <ss.NS; ++s) {
    for (int a = 0; a < pars.base.options.p_idx_fertility_first; ++a) {
      for (int hd = 0; hd < hc_ss.hc1DS; ++hd) {
        for (int cat = 0; cat < hc_ss.hcTT; ++cat) {
          if (cpars.hc_art_val(time_step) > 0) {
            intermediate.children.hc_art_init(hd, cat, a, s) += intermediate.children.hc_art_need(hd, cat, a, s);
            for (int dur = 0; dur < ss.hTS; ++dur) {
              if (intermediate.children.hc_art_init(hd, cat, a, s) < 0) {
                intermediate.children.hc_art_init(hd, cat, a, s) = 0.0;
              }else{
                intermediate.children.hc_art_init(hd, cat, a, s) = intermediate.children.hc_art_init(hd, cat, a, s);
              }
            }// end ss.hTS
          }
        }// end hcTT
      }// end hc_ss.hc1DS
    }// end a
  }// end ss.NS

  for (int s = 0; s <ss.NS; ++s) {
    for (int a = 0; a < pars.base.options.p_idx_fertility_first; ++a) {
      for (int hd = 0; hd < hc_ss.hc1DS; ++hd) {
        for (int cat = 0; cat < hc_ss.hcTT; ++cat) {
          intermediate.children.hc_art_init_total += intermediate.children.hc_art_init(hd, cat, a, s);
        }// end hcTT
      }// end hc_ss.hc1DS
    }// end a
  }// end ss.NS

  //!!! TODO: fix order of for loop
  for (int s = 0; s <ss.NS; ++s) {
    for (int hd = 0; hd < hc_ss.hc1DS; ++hd) {
      for (int a = 0; a < pars.base.options.p_idx_fertility_first; ++a) {
        intermediate.children.hc_death_rate = 0.0;
        intermediate.children.hc_art_grad(0, hd, a, s) = 0.0;
        if (a < hc_ss.hc2_agestart) {
          if (state_next.children.hc1_art_pop(0, hd, a, s) >0) {
            intermediate.children.hc_death_rate =  cpars.hc_art_mort_rr(0, a, time_step) * 0.5 * (cpars.hc1_art_mort(hd, 0, a) + cpars.hc1_art_mort(hd, 1, a));
            state_next.children.hc1_art_aids_deaths(0,hd, a, s) =  intermediate.children.hc_death_rate * state_next.children.hc1_art_pop(0, hd, a, s);
            intermediate.children.hc_art_grad(0,hd, a, s) -= state_next.children.hc1_art_aids_deaths(0,hd, a, s);
            state_next.children.hc1_art_pop(0, hd,  a, s) += intermediate.children.hc_art_grad(0, hd, a, s);
          }
        } else if (hd < hc_ss.hc2DS && state_next.children.hc2_art_pop(0, hd, a-hc_ss.hc2_agestart, s) >0) {
          intermediate.children.hc_death_rate =  cpars.hc_art_mort_rr(0, a, time_step) * 0.5 * (cpars.hc2_art_mort(hd, 0, a-hc_ss.hc2_agestart) + cpars.hc2_art_mort(hd, 1, a-hc_ss.hc2_agestart));
          state_next.children.hc2_art_aids_deaths(0,hd, a-hc_ss.hc2_agestart, s) =  intermediate.children.hc_death_rate * state_next.children.hc2_art_pop(0, hd, a-hc_ss.hc2_agestart, s);
          intermediate.children.hc_art_grad(0, hd, a, s) -= state_next.children.hc2_art_aids_deaths(0,hd, a-hc_ss.hc2_agestart, s);
          state_next.children.hc2_art_pop(0, hd,  a-hc_ss.hc2_agestart, s) += intermediate.children.hc_art_grad(0, hd, a, s);
        }
      }// end a
    }// end hc_ss.hc1DS
  }// end ss.NS

  for (int s = 0; s <ss.NS; ++s) {
    for (int hd = 0; hd < hc_ss.hc1DS; ++hd) {
      for (int a = 0; a < pars.base.options.p_idx_fertility_first; ++a) {
        intermediate.children.hc_death_rate = 0;
        intermediate.children.hc_art_grad(2, hd, a, s) = 0.0;

        if (a < hc_ss.hc2_agestart) {
          intermediate.children.hc_death_rate = cpars.hc_art_mort_rr(2, a, time_step) * cpars.hc1_art_mort(hd, 2, a);
          state_next.children.hc1_art_aids_deaths(2,hd, a, s) =  state_next.children.hc1_art_pop(2, hd, a, s) * intermediate.children.hc_death_rate;
          intermediate.children.hc_art_grad(2,hd, a, s) -= state_next.children.hc1_art_aids_deaths(2,hd, a, s);
          state_next.children.hc1_art_pop(2, hd,  a, s) += intermediate.children.hc_art_grad(2, hd, a, s);
        } else if (hd < (hc_ss.hc2DS)) {
          intermediate.children.hc_death_rate = cpars.hc_art_mort_rr(2, a, time_step) * cpars.hc2_art_mort(hd, 2, a-hc_ss.hc2_agestart);
          state_next.children.hc2_art_aids_deaths(2,hd, a-hc_ss.hc2_agestart, s) =  state_next.children.hc2_art_pop(2, hd, a-hc_ss.hc2_agestart, s) * intermediate.children.hc_death_rate;
          intermediate.children.hc_art_grad(2,hd, a, s) -= state_next.children.hc2_art_aids_deaths(2,hd, a-hc_ss.hc2_agestart, s);
          state_next.children.hc2_art_pop(2, hd,  a-hc_ss.hc2_agestart, s) += intermediate.children.hc_art_grad(2, hd, a, s);
        }
      }// end a
    }// end hc_ss.hc1DS
  }// end ss.NS
}

template<typename ModelVariant, typename real_type>
void hc_art_num_num(int time_step,
                    const Parameters<ModelVariant, real_type> &pars,
                    const State<ModelVariant, real_type> &state_curr,
                    State<ModelVariant, real_type> &state_next,
                    IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "run_hiv_child_infections can only be called for model variants where run_child_model is true");
  constexpr auto ss = StateSpace<ModelVariant>().base;
  constexpr auto hc_ss = StateSpace<ModelVariant>().children;
  const auto cpars = pars.children.children;

  //Remove how many that are already on ART
  state_next.children.hc_art_num =  (cpars.hc_art_val(time_step) + cpars.hc_art_val(time_step-1)) / 2 ;
  for (int s = 0; s <ss.NS; ++s) {
    for (int a = 0; a < pars.base.options.p_idx_fertility_first; ++a) {
      for (int hd = 0; hd < hc_ss.hc1DS; ++hd) {
        for (int dur = 0; dur < ss.hTS; ++dur) {
          if (a < hc_ss.hc2_agestart) {
            state_next.children.hc_art_num -= state_next.children.hc1_art_pop(dur, hd, a, s);
          }else if (hd < (hc_ss.hc2DS)) {
            state_next.children.hc_art_num -= state_next.children.hc2_art_pop(dur, hd, a-hc_ss.hc2_agestart, s);
          }
        }// end ss.hTS
      }// end hc_ss.hc1DS
    }// end a
  }// end ss.NS

  if (intermediate.children.hc_art_init_total < state_next.children.hc_art_num) {
    state_next.children.hc_art_num = intermediate.children.hc_art_init_total;
  }
  if (state_next.children.hc_art_num < 0) {
    state_next.children.hc_art_num =  0;
  }
}

template<typename ModelVariant, typename real_type>
void hc_art_pct_pct(int time_step,
                    const Parameters<ModelVariant, real_type> &pars,
                    const State<ModelVariant, real_type> &state_curr,
                    State<ModelVariant, real_type> &state_next,
                    IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "run_hiv_child_infections can only be called for model variants where run_child_model is true");
  constexpr auto ss = StateSpace<ModelVariant>().base;
  constexpr auto hc_ss = StateSpace<ModelVariant>().children;
  const auto cpars = pars.children.children;

  for (int s = 0; s <ss.NS; ++s) {
    for (int a = 0; a < pars.base.options.p_idx_fertility_first; ++a) {
      for (int hd = 0; hd < hc_ss.hc1DS; ++hd) {
        for (int dur = 0; dur < ss.hTS; ++dur) {
          if (a < hc_ss.hc2_agestart) {
            state_next.children.hc_art_num += state_next.children.hc1_art_pop(dur, hd, a, s);
            state_next.children.hc_art_num += state_next.children.hc1_art_aids_deaths(dur,hd, a, s);
          }else if (hd < (hc_ss.hc2DS )) {
            state_next.children.hc_art_num += state_next.children.hc2_art_pop(dur, hd, a-hc_ss.hc2_agestart, s);
            state_next.children.hc_art_num += state_next.children.hc2_art_aids_deaths(dur,hd, a-hc_ss.hc2_agestart, s);
          }
        } //end ss.hTS
      } // end ss.hC1_disease_stages
    } // end a
  } // end ss.NS
  state_next.children.hc_art_num =  state_next.children.hc_art_num * (cpars.hc_art_val(time_step) + cpars.hc_art_val(time_step-1)) / 2;

  //Remove how many that are already on ART
  for (int s = 0; s <ss.NS; ++s) {
    for (int a = 0; a < pars.base.options.p_idx_fertility_first; ++a) {
      for (int hd = 0; hd < hc_ss.hc1DS; ++hd) {
        for (int dur = 0; dur < ss.hTS; ++dur) {
          if (a < hc_ss.hc2_agestart) {
            state_next.children.hc_art_num -= state_next.children.hc1_art_pop(dur, hd, a, s);
          }else if (hd < (hc_ss.hc2DS)) {
            state_next.children.hc_art_num -= state_next.children.hc2_art_pop(dur, hd, a-hc_ss.hc2_agestart, s);
          }
        } // end ss.hTS
      } // end ss.hC1_disease_stages
    } // end a
  } // end ss.NS
  if (intermediate.children.hc_art_init_total < state_next.children.hc_art_num) {
    state_next.children.hc_art_num = intermediate.children.hc_art_init_total;
  }
}

template<typename ModelVariant, typename real_type>
void hc_art_num_pct(int time_step,
                    const Parameters<ModelVariant, real_type> &pars,
                    const State<ModelVariant, real_type> &state_curr,
                    State<ModelVariant, real_type> &state_next,
                    IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "run_hiv_child_infections can only be called for model variants where run_child_model is true");
  constexpr auto ss = StateSpace<ModelVariant>().base;
  constexpr auto hc_ss = StateSpace<ModelVariant>().children;
  const auto cpars = pars.children.children;

  //Remove how many that are already on ART
  for (int s = 0; s <ss.NS; ++s) {
    for (int a = 0; a < pars.base.options.p_idx_fertility_first; ++a) {
      for (int hd = 0; hd < hc_ss.hc1DS; ++hd) {
        for (int dur = 0; dur < ss.hTS; ++dur) {
          if (a < hc_ss.hc2_agestart) {
            state_next.children.hc_art_num += state_next.children.hc1_art_pop(dur, hd, a, s) +  state_next.children.hc1_art_aids_deaths(dur,hd, a, s);
          } else if (hd < (hc_ss.hc2DS)) {
            state_next.children.hc_art_num += state_next.children.hc2_art_pop(dur, hd, a-hc_ss.hc2_agestart, s) + state_next.children.hc2_art_aids_deaths(dur,hd, a-hc_ss.hc2_agestart, s);
          }
        } // end ss.hTS
      } // end hc_ss.hc1DS
    } //end a
  } //end ss.NS
  state_next.children.hc_art_num = (cpars.hc_art_val(time_step-1) + (state_next.children.hc_art_num * cpars.hc_art_val(time_step))) / 2;

  for (int s = 0; s <ss.NS; ++s) {
    for (int a = 0; a < pars.base.options.p_idx_fertility_first; ++a) {
      for (int hd = 0; hd < hc_ss.hc1DS; ++hd) {
        for (int dur = 0; dur < ss.hTS; ++dur) {
          if (a < hc_ss.hc2_agestart) {
            state_next.children.hc_art_num -= state_next.children.hc1_art_pop(dur, hd, a, s);
          } else if (hd < (hc_ss.hc2DS )) {
            state_next.children.hc_art_num -= state_next.children.hc2_art_pop(dur, hd, a-hc_ss.hc2_agestart, s);
          }
        } // end ss.hTS
      } // end hc_ss.hc1DS
    } // end a
  } // end ss.NS

  if (state_next.children.hc_art_num < 0) {
    state_next.children.hc_art_num = 0;
  }
  if (intermediate.children.hc_art_init_total < state_next.children.hc_art_num) {
    state_next.children.hc_art_num = intermediate.children.hc_art_init_total;
  }
}

template<typename ModelVariant, typename real_type>
void hc_art_pct_num(int time_step,
                    const Parameters<ModelVariant, real_type> &pars,
                    const State<ModelVariant, real_type> &state_curr,
                    State<ModelVariant, real_type> &state_next,
                    IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "run_hiv_child_infections can only be called for model variants where run_child_model is true");
  constexpr auto ss = StateSpace<ModelVariant>().base;
  constexpr auto hc_ss = StateSpace<ModelVariant>().children;
  const auto cpars = pars.children.children;

  for (int s = 0; s <ss.NS; ++s) {
    for (int a = 0; a < pars.base.options.p_idx_fertility_first; ++a) {
      for (int hd = 0; hd < hc_ss.hc1DS; ++hd) {
        for (int dur = 0; dur < ss.hTS; ++dur) {
          if (a < hc_ss.hc2_agestart) {
            state_next.children.hc_art_num -= state_next.children.hc1_art_pop(dur, hd, a, s);
          }else if (hd < (hc_ss.hc2DS)) {
            state_next.children.hc_art_num -= state_next.children.hc2_art_pop(dur, hd, a-hc_ss.hc2_agestart, s);
          }
        } // end ss.hTS
      } //end hc_ss.hc1DS
    } // end a
  } //end ss.NS

  state_next.children.hc_art_num = (state_curr.children.hc_art_num + cpars.hc_art_val(time_step)) / 2;

  //Remove how many that are already on ART
  for (int s = 0; s <ss.NS; ++s) {
    for (int a = 0; a < pars.base.options.p_idx_fertility_first; ++a) {
      for (int hd = 0; hd < hc_ss.hc1DS; ++hd) {
        for (int dur = 0; dur < ss.hTS; ++dur) {
          if (a < hc_ss.hc2_agestart) {
            state_next.children.hc_art_num -= state_next.children.hc1_art_pop(dur, hd, a, s);
          } else if (hd < (hc_ss.hc2DS)) {
            state_next.children.hc_art_num -= state_next.children.hc2_art_pop(dur, hd, a-hc_ss.hc2_agestart, s);
          }
        } //end ss.hTS
      } //end hc_ss.hc1DS
    } //end a
  } //end ss.NS

  if (state_next.children.hc_art_num < 0) {
    state_next.children.hc_art_num = 0;
  }
  if (intermediate.children.hc_art_init_total < state_next.children.hc_art_num) {
    state_next.children.hc_art_num = intermediate.children.hc_art_init_total;
  }
}


template<typename ModelVariant, typename real_type>
void hc_art_initiation_by_five_year_age_groups(int time_step,
                              const Parameters<ModelVariant, real_type> &pars,
                              const State<ModelVariant, real_type> &state_curr,
                              State<ModelVariant, real_type> &state_next,
                              IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "run_hiv_child_infections can only be called for model variants where run_child_model is true");
  constexpr auto ss = StateSpace<ModelVariant>().base;
  constexpr auto hc_ss = StateSpace<ModelVariant>().children;
  const auto cpars = pars.children.children;

  // for (int s = 0; s <ss.NS; ++s) {
  //   for (int a = 0; a < pars.base.options.p_idx_fertility_first; ++a) {
  //     for (int hd = 0; hd < hc_ss.hc1DS; ++hd) {
  //       for (int cat = 0; cat < hc_ss.hcTT; ++cat) {
  //         intermediate.children.hc_initByAge += intermediate.children.hc_art_init(hd, cat, a, s) * cpars.hc_art_init_dist(a, time_step);
  //       } //end hcTT
  //     } // end hc_ss.hc1DS
  //   } //end a
  // } // end ss.NS
  // if (intermediate.children.hc_initByAge == 0.0) {
  //   intermediate.children.hc_adj = 1.0 ;
  // }else{
  //   intermediate.children.hc_adj = state_next.children.hc_art_num / intermediate.children.hc_initByAge;
  // }
  // for (int s = 0; s <ss.NS; ++s) {
  //   for (int cat = 0; cat < hc_ss.hcTT; ++cat) {
  //     for (int a = 0; a < pars.base.options.p_idx_fertility_first; ++a) {
  //       for (int hd = 0; hd < hc_ss.hc1DS; ++hd) {
  //         if ((intermediate.children.hc_adj * cpars.hc_art_init_dist(a, time_step)) > 1.0) {
  //           intermediate.children.hc_art_scalar = 1.0;
  //         }else{
  //           intermediate.children.hc_art_scalar = intermediate.children.hc_adj * cpars.hc_art_init_dist(a, time_step);
  //         }
  //         if (state_next.children.hc_art_num > 0.0) {
  //           intermediate.children.hc_art_scalar = intermediate.children.hc_art_scalar;
  //         }else{
  //           intermediate.children.hc_art_scalar = 0.0;
  //         }
  //         if (a < hc_ss.hc2_agestart) {
  //           state_next.children.hc1_art_pop(0, hd, a, s) += intermediate.children.hc_art_scalar * intermediate.children.hc_art_init(hd, cat, a, s);
  //         }else if (hd < (hc_ss.hc2DS)) {
  //           state_next.children.hc2_art_pop(0, hd, a - hc_ss.hc2_agestart, s) += intermediate.children.hc_art_scalar * intermediate.children.hc_art_init(hd, cat, a, s);
  //         }
  //         if (a < hc_ss.hc2_agestart) {
  //           state_next.children.hc1_hiv_pop(hd, cat, a, s) -= intermediate.children.hc_art_scalar * intermediate.children.hc_art_init(hd, cat, a, s);
  //         }else if (hd < (hc_ss.hc2DS )) {
  //           state_next.children.hc2_hiv_pop(hd, cat, a - hc_ss.hc2_agestart, s) -=  intermediate.children.hc_art_scalar * intermediate.children.hc_art_init(hd, cat, a, s);
  //         }
  //
  //       } //end hc_ss.hc1DS
  //     } // end a
  //   } // end hcTT
  // } // end ss.NS
}


template<typename ModelVariant, typename real_type>
void run_child_art_initiation(int time_step,
                              const Parameters<ModelVariant, real_type> &pars,
                              const State<ModelVariant, real_type> &state_curr,
                              State<ModelVariant, real_type> &state_next,
                              IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "run_hiv_child_infections can only be called for model variants where run_child_model is true");
  constexpr auto ss = StateSpace<ModelVariant>().base;
  constexpr auto hc_ss = StateSpace<ModelVariant>().children;
  const auto cpars = pars.children.children;


  internal::hc_initiate_art_by_age(time_step, pars, state_curr, state_next, intermediate); //can be used from
  internal::hc_initiate_art_by_cd4(time_step, pars, state_curr, state_next, intermediate);
  internal::hc_adjust_art_initiates_for_mort_age_spec(time_step, pars, state_curr, state_next, intermediate);

  //Progress ART to the correct time on ART
  for (int hd = 0; hd < hc_ss.hc1DS; ++hd) {
    for (int a = 0; a < pars.base.options.p_idx_fertility_first; ++a) {
      for (int s = 0; s <ss.NS; ++s) {
        if (a < hc_ss.hc2_agestart) {
          if (state_next.children.hc1_art_pop(0, hd, a, s) > 0) {
            state_next.children.hc1_art_pop(1, hd, a, s) += state_next.children.hc1_art_pop(0, hd, a, s);
            state_next.children.hc1_art_pop(0, hd, a, s) -= state_next.children.hc1_art_pop(0, hd, a, s);
          }
        }else if (hd < (hc_ss.hc2DS)) {
            state_next.children.hc2_art_pop(1, hd, a-hc_ss.hc2_agestart, s) += state_next.children.hc2_art_pop(0, hd, a-hc_ss.hc2_agestart, s);
            state_next.children.hc2_art_pop(0, hd, a-hc_ss.hc2_agestart, s) -= state_next.children.hc2_art_pop(0, hd, a-hc_ss.hc2_agestart, s);
          }
      }//end ss.NS
    }// end a
  }// end hc_ss.hc1DS

  if (!cpars.hc_art_isperc(time_step) && !cpars.hc_art_isperc(time_step-1)) { // both numbers
   internal::hc_art_num_num(time_step, pars, state_curr, state_next, intermediate);
  } else if (cpars.hc_art_isperc(time_step) && cpars.hc_art_isperc(time_step-1)) { // both percentages
    internal::hc_art_pct_pct(time_step, pars, state_curr, state_next, intermediate);
  } else if (cpars.hc_art_isperc(time_step) && !cpars.hc_art_isperc(time_step-1)) { // num to percentage
    internal::hc_art_num_pct(time_step, pars, state_curr, state_next, intermediate);
  } else if (cpars.hc_art_isperc(time_step-1) && !cpars.hc_art_isperc(time_step)) { //percentage to num
    internal::hc_art_pct_num(time_step, pars, state_curr, state_next, intermediate);
  }

  internal::hc_art_initiation_by_age(time_step, pars, state_curr, state_next, intermediate);
}



}// namespace internal

template<typename ModelVariant, typename real_type>
void run_art_initiation_by_five_year_age_group(int time_step,
                                const Parameters<ModelVariant, real_type> &pars,
                                const State<ModelVariant, real_type> &state_curr,
                                State<ModelVariant, real_type> &state_next,
                                internal::IntermediateData<ModelVariant, real_type> &intermediate) {

  internal::run_child_art_initiation(time_step, pars, state_curr, state_next, intermediate);
  internal::run_child_art_mortality(time_step, pars, state_curr, state_next, intermediate);
  }

} // namespace leapfrog
