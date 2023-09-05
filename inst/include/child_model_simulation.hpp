#pragma once

#include "types.hpp"

namespace leapfrog {

namespace internal {

template<typename ModelVariant, typename real_type>
void run_child_ageing(int time_step,
                      const Parameters<ModelVariant, real_type> &pars,
                      const State<ModelVariant, real_type> &state_curr,
                      State<ModelVariant, real_type> &state_next,
                      IntermediateData<ModelVariant, real_type> &intermediate) {
  const auto demog = pars.base.demography;
  const auto cpars = pars.children.children;
  static_assert(ModelVariant::run_child_model,
                "run_hiv_child_infections can only be called for model variants where run_child_model is true");
  constexpr auto ss = StateSpace<ModelVariant>().base;
  constexpr auto hc_ss = StateSpace<ModelVariant>().children;

  for (int s = 0; s < ss.NS; ++s) {
    //less than 5 because there is a cd4 transition between ages 4 and 5

    for (int a = 1; a < hc_ss.hc2_agestart; ++a) {
      for (int hd = 0; hd < hc_ss.hc1DS; ++hd) {
         for (int cat = 0; cat < hc_ss.hcTT; ++cat) {
          state_next.children.hc1_hiv_pop(hd, cat, a, s) +=
              state_curr.children.hc1_hiv_pop(hd, cat, a - 1, s) * demog.survival_probability(a, s, time_step);
        }
          for (int dur = 0; dur < ss.hTS; ++dur) {
          state_next.children.hc1_art_pop(dur, hd, a, s) +=
              state_curr.children.hc1_art_pop(dur, hd, a - 1, s) * demog.survival_probability(a, s, time_step);
        }
      }
    }
  }

  for (int s = 0; s < ss.NS; ++s) {

    for (int hd = 0; hd < hc_ss.hc1DS; ++hd) {
      for (int hd_alt = 0; hd_alt < hc_ss.hc2DS; ++hd_alt) {
        for (int cat = 0; cat < hc_ss.hcTT; ++cat) {
          state_next.children.hc2_hiv_pop(hd_alt, cat, 0, s) += state_curr.children.hc1_hiv_pop(hd, cat, hc_ss.hc1_ageend, s) * demog.survival_probability(hc_ss.hc2_agestart, s, time_step) * cpars.hc_cd4_transition(hd_alt, hd);
        }
          for (int dur = 0; dur < ss.hTS; ++dur) {
          state_next.children.hc2_art_pop(dur, hd_alt, 0, s) += state_curr.children.hc1_art_pop(dur, hd, hc_ss.hc1_ageend, s) * demog.survival_probability(hc_ss.hc2_agestart, s, time_step) * cpars.hc_cd4_transition(hd_alt, hd);
        }
      }
    }
  }

  for (int s = 0; s < ss.NS; ++s) {
    for (int a = (hc_ss.hc2_agestart + 1); a < pars.base.options.p_idx_fertility_first; ++a) {
      for (int hd = 0; hd < hc_ss.hc2DS; ++hd) {
       for (int cat = 0; cat < hc_ss.hcTT; ++cat) {
          state_next.children.hc2_hiv_pop(hd, cat, a - hc_ss.hc2_agestart, s) +=
              state_curr.children.hc2_hiv_pop(hd, cat, a - hc_ss.hc2_agestart - 1, s) *
              demog.survival_probability(a, s, time_step);
        }
          for (int dur = 0; dur < ss.hTS; ++dur) {
          state_next.children.hc2_art_pop(dur, hd, a - hc_ss.hc2_agestart, s) +=
              state_curr.children.hc2_art_pop(dur, hd, a - hc_ss.hc2_agestart - 1, s) *
              demog.survival_probability(a, s, time_step);
        }
      }
    }
  }


}

template<typename ModelVariant, typename real_type>
void run_child_hiv_infections(int time_step,
                              const Parameters<ModelVariant, real_type> &pars,
                              const State<ModelVariant, real_type> &state_curr,
                              State<ModelVariant, real_type> &state_next,
                              IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "run_hiv_child_infections can only be called for model variants where run_child_model is true");
  constexpr auto ss = StateSpace<ModelVariant>().base;
  constexpr auto hc_ss = StateSpace<ModelVariant>().children;
  const auto cpars = pars.children.children;

  for (int s = 0; s < ss.NS; ++s) {
    // Run only first 5 age groups in total population 0, 1, 2, 3, 4
    for (int a = 0; a < hc_ss.hc2_agestart; ++a) {
      if (cpars.hc_nosocomial(time_step) > 0) {
        //5.0 is used because we want to evenly distribute across the 5 age groups in 0-4
        state_next.base.p_infections(a, s) = cpars.hc_nosocomial(time_step) / (5.0 * ss.NS);
        state_next.base.p_hiv_pop(a, s) += state_next.base.p_infections(a, s);

        for (int hd = 0; hd < hc_ss.hc1DS; ++hd) {
          // putting them all in perinatal hcTT to match spec nosocomial
          if (cpars.hc1_cd4_dist(hd) > 0) {
            state_next.children.hc1_hiv_pop(hd, 0, a, s) += state_next.base.p_infections(a, s) * cpars.hc1_cd4_dist(hd);

          }
        }
      }
    }
  }
}


template<typename ModelVariant, typename real_type>
void run_child_natural_history(int time_step,
                               const Parameters<ModelVariant, real_type> &pars,
                               const State<ModelVariant, real_type> &state_curr,
                              State<ModelVariant, real_type> &state_next,
                              IntermediateData<ModelVariant, real_type> &intermediate) {
 static_assert(ModelVariant::run_child_model,
               "run_hiv_child_infections can only be called for model variants where run_child_model is true");
 constexpr auto ss = StateSpace<ModelVariant>().base;
 constexpr auto hc_ss = StateSpace<ModelVariant>().children;
 const auto cpars = pars.children.children;


  for (int s = 0; s < ss.NS; ++s) {
    for (int a = 0; a < hc_ss.hc2_agestart; ++a) {
      for (int cat = 0; cat < hc_ss.hcTT; ++cat) {
          for (int hd = 0; hd < hc_ss.hc1DS; ++hd) {
          intermediate.children.hc_posthivmort(hd, cat, a, s) = state_next.children.hc1_hiv_pop(hd, cat, a, s) -
            (1.0 - cpars.ctx_effect * cpars.ctx_val(time_step)) * state_next.children.hc1_hiv_pop(hd, cat, a, s) * cpars.hc1_cd4_mort(hd, cat, a);

        }
      }
    }
  }

  for (int s = 0; s < ss.NS; ++s) {
    for (int a = hc_ss.hc2_agestart; a < pars.base.options.p_idx_fertility_first; ++a) {
      for (int cat = 0; cat < hc_ss.hcTT; ++cat) {
        for (int hd = 0; hd < hc_ss.hc2DS; ++hd) {
          intermediate.children.hc_posthivmort(hd, cat, a, s) = state_next.children.hc2_hiv_pop(hd, cat, a - hc_ss.hc2_agestart, s) -
            (1.0 - cpars.ctx_effect * cpars.ctx_val(time_step)) * state_next.children.hc2_hiv_pop(hd, cat, a - hc_ss.hc2_agestart, s) * cpars.hc2_cd4_mort(hd, cat, a - hc_ss.hc2_agestart);
        }
      }
    }
  }

  //progress through CD4 categories
  for (int s = 0; s < ss.NS; ++s) {
    for (int a = 0; a < hc_ss.hc2_agestart; ++a) {
      for (int cat = 0; cat < hc_ss.hcTT; ++cat) {
        for (int hd = 1; hd < hc_ss.hc1DS; ++hd) {
          intermediate.children.hc_grad(hd - 1, cat, a, s) -=
              (intermediate.children.hc_posthivmort(hd - 1, cat, a, s) * cpars.hc1_cd4_prog(hd - 1) +
               state_next.children.hc1_hiv_pop(hd - 1, cat, a, s) * cpars.hc1_cd4_prog(hd - 1)) /
              2; //moving to next cd4 category
          intermediate.children.hc_grad(hd, cat, a, s) +=
              (intermediate.children.hc_posthivmort(hd - 1, cat, a, s) * cpars.hc1_cd4_prog(hd - 1) +
               state_next.children.hc1_hiv_pop(hd - 1, cat, a, s) * cpars.hc1_cd4_prog(hd - 1)) /
              2; //moving into this cd4 category


 for (int s = 0; s < ss.NS; ++s) {
   for (int hd = 0; hd < hc_ss.hc2DS; ++hd) {
     for (int a = hc_ss.hc2_agestart; a < pars.base.options.p_idx_fertility_first; ++a) {
       for (int cat = 0; cat < hc_ss.hcTT; ++cat) {
         intermediate.children.hc_posthivmort(hd, cat, a, s) += state_next.children.hc2_hiv_pop(hd, cat, a - hc_ss.hc2_agestart, s) - (1 - cpars.ctx_effect * cpars.ctx_val(time_step)) * state_curr.children.hc2_hiv_pop(hd, cat, a - hc_ss.hc2_agestart, s) * cpars.hc2_cd4_mort(hd, cat, a - hc_ss.hc2_agestart);
         state_next.children.hc2_noart_aids_deaths(hd, cat, a, s) +=  (1 - cpars.ctx_effect * cpars.ctx_val(time_step)) * state_curr.children.hc2_hiv_pop(hd, cat, a, s) * cpars.hc2_cd4_mort(hd, cat, a - 5); // output hiv deaths, aggregated across transmission category
       }
     }
   }
 }

 //progress through CD4 categories
 for (int s = 0; s < ss.NS; ++s) {
   for (int hd = 1; hd < hc_ss.hc1DS; ++hd) {
     for (int a = 0; a < hc_ss.hc2_agestart; ++a) {
       for (int cat = 0; cat < hc_ss.hcTT; ++cat) {
         intermediate.children.hc_grad(hd - 1, cat, a, s) -=  (intermediate.children.hc_posthivmort(hd - 1, cat, a, s) * cpars.hc1_cd4_prog(hd - 1) + state_next.children.hc1_hiv_pop(hd - 1, cat, a, s) * cpars.hc1_cd4_prog(hd - 1)) / 2; //moving to next cd4 category
         intermediate.children.hc_grad(hd, cat, a, s) += (intermediate.children.hc_posthivmort(hd - 1, cat, a, s) * cpars.hc1_cd4_prog(hd - 1) + state_next.children.hc1_hiv_pop(hd - 1, cat, a, s) * cpars.hc1_cd4_prog(hd - 1)) / 2; //moving into this cd4 category
       }
     }
   }
 }

 //progress through CD4 categories
 for (int s = 0; s < ss.NS; ++s) {
   for (int hd = 1; hd < hc_ss.hc2DS; ++hd) {
     for (int a = hc_ss.hc2_agestart; a < pars.base.options.p_idx_fertility_first; ++a) {
       for (int cat = 0; cat < hc_ss.hcTT; ++cat) {
         intermediate.children.hc_grad(hd - 1, cat, a, s) -= (intermediate.children.hc_posthivmort(hd - 1, cat, a, s) * cpars.hc2_cd4_prog(hd - 1) + state_next.children.hc2_hiv_pop(hd - 1, cat, a - hc_ss.hc2_agestart, s) * cpars.hc2_cd4_prog(hd - 1)) / 2; //moving to next cd4 category
         intermediate.children.hc_grad(hd, cat, a, s) += (intermediate.children.hc_posthivmort(hd - 1, cat, a, s) * cpars.hc2_cd4_prog(hd - 1) + state_next.children.hc2_hiv_pop(hd - 1, cat, a - hc_ss.hc2_agestart, s) * cpars.hc2_cd4_prog(hd - 1)) / 2; //moving into this cd4 category
       }
     }
   }
 }
}

template<typename ModelVariant, typename real_type>
void run_child_hiv_mort(int time_step,
                        const Parameters<ModelVariant, real_type> &pars,
                               const State<ModelVariant, real_type> &state_curr,
                               State<ModelVariant, real_type> &state_next,
                               IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "run_hiv_child_infections can only be called for model variants where run_child_model is true");
  constexpr auto ss = StateSpace<ModelVariant>().base;
  constexpr auto hc_ss = StateSpace<ModelVariant>().children;
  const auto cpars = pars.children.children;

  for (int s = 0; s < ss.NS; ++s) {
    for (int a = 0; a < hc_ss.hc2_agestart; ++a) {
      for (int cat = 0; cat < hc_ss.hcTT; ++cat) {
          for (int hd = 0; hd < hc_ss.hc1DS; ++hd) {
            intermediate.children.hc_grad(hd, cat, a, s) -=
              state_next.children.hc1_hiv_pop(hd, cat, a, s) * cpars.hc1_cd4_mort(hd, cat, a);

         }
      }
    }
  }

  for (int s = 0; s < ss.NS; ++s) {
    for (int a = hc_ss.hc2_agestart; a < pars.base.options.p_idx_fertility_first; ++a) {
      for (int cat = 0; cat < hc_ss.hcTT; ++cat) {
            for (int hd = 0; hd < hc_ss.hc2DS; ++hd) {
            intermediate.children.hc_grad(hd, cat, a, s) -=
              state_next.children.hc2_hiv_pop(hd, cat, a - hc_ss.hc2_agestart, s) *
              cpars.hc2_cd4_mort(hd, cat, a - hc_ss.hc2_agestart);
          }
      }
    }
  }
}

template<typename ModelVariant, typename real_type>
void add_child_grad(int time_step,
                    const Parameters<ModelVariant, real_type> &pars,
                    const State<ModelVariant, real_type> &state_curr,
                    State<ModelVariant, real_type> &state_next,
                    IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "run_hiv_child_infections can only be called for model variants where run_child_model is true");
  constexpr auto ss = StateSpace<ModelVariant>().base;
  constexpr auto hc_ss = StateSpace<ModelVariant>().children;

    //add on transitions
    for (int s = 0; s < ss.NS; ++s) {
      for (int a = 0; a < hc_ss.hc2_agestart; ++a) {
        for (int cat = 0; cat < hc_ss.hcTT; ++cat) {
            for (int hd = 0; hd < hc_ss.hc1DS; ++hd) {
              state_next.children.hc1_hiv_pop(hd, cat, a, s) += intermediate.children.hc_grad(hd, cat, a, s);
            }// end hc1DS
        }//end cat
        }//end a
      }// end s

      for (int s = 0; s < ss.NS; ++s) {
        for (int a = hc_ss.hc2_agestart; a < pars.base.options.p_idx_fertility_first; ++a) {
          for (int cat = 0; cat < hc_ss.hcTT; ++cat) {
              for (int hd = 0; hd < hc_ss.hc2DS; ++hd) {
                state_next.children.hc2_hiv_pop(hd, cat, a - hc_ss.hc2_agestart, s) += intermediate.children.hc_grad(hd, cat, a, s);
              }//end hc2DS
          }//end cat
          }//end a
        }// end s

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


  //all children under a certain age eligible for ART
  for(int s = 0; s <ss.NS; ++s)  {
    for(int cat = 0; cat < hc_ss.hcTT; ++cat) {
      for(int a = 0; a < pars.base.options.p_idx_fertility_first; ++a) {
        for(int hd = 0; hd < hc_ss.hc1DS; ++hd) {
          if(a < cpars.hc_art_elig_age(time_step)){
            if(a < hc_ss.hc2_agestart){
              intermediate.children.hc_art_need(hd, cat, a, s) += state_next.children.hc1_hiv_pop(hd, cat, a, s);
            }else{
            if(hd < hc_ss.hc2DS){
              intermediate.children.hc_art_need(hd, cat, a, s) += state_next.children.hc2_hiv_pop(hd, cat, a - hc_ss.hc2_agestart, s);
            }
            }
          }
        } // end ss.hc1DS
      } // end a
    } // end hcTT
  } // end ss.NS

  //all children under a certain CD4 eligible for ART
  for(int s = 0; s <ss.NS; ++s)  {
    for(int cat = 0; cat < hc_ss.hcTT; ++cat) {
      for(int a = cpars.hc_art_elig_age(time_step); a < pars.base.options.p_idx_fertility_first; ++a) {
        for(int hd = 0; hd < hc_ss.hc1DS; ++hd) {
          if(hd > (cpars.hc_art_elig_cd4(a, time_step) - 2)){ //-2 to account for the zero-based indexing and basically as >=
            if(a < hc_ss.hc2_agestart){
              intermediate.children.hc_art_need(hd, cat, a, s) += state_next.children.hc1_hiv_pop(hd, cat, a, s);
            }else{
              if(hd < hc_ss.hc2DS){
                intermediate.children.hc_art_need(hd, cat, a, s) += state_next.children.hc2_hiv_pop(hd, cat, a - hc_ss.hc2_agestart, s);
              }
            }
          }
        } // end ss.hc1DS
      } // end a
    } // end hcTT
  } // end ss.NS

  for(int s = 0; s <ss.NS; ++s) {
    for(int a = 0; a < pars.base.options.p_idx_fertility_first; ++a) {
      for(int hd = 0; hd < hc_ss.hc1DS; ++hd) {
        for(int cat = 0; cat < hc_ss.hcTT; ++cat) {
          state_next.children.hc_art_num += intermediate.children.hc_art_need(hd, cat, a, s);
        } // end hcTT
      } // end ss.hc1DS
    } // end a
  } // end ss.NS

  //how many should initialize ART
  for(int s = 0; s <ss.NS; ++s) {
    for(int a = 0; a < pars.base.options.p_idx_fertility_first; ++a) {
      for(int hd = 0; hd < hc_ss.hc1DS; ++hd) {
        for(int cat = 0; cat < hc_ss.hcTT; ++cat) {
          if(cpars.hc_art_val(time_step) > 0){
            intermediate.children.hc_art_init(hd, cat, a, s) += intermediate.children.hc_art_need(hd, cat, a, s);
            for(int dur = 0; dur < ss.hTS; ++dur) {
              if(intermediate.children.hc_art_init(hd, cat, a, s) < 0){
                intermediate.children.hc_art_init(hd, cat, a, s) = 0.0;
              }else{
                intermediate.children.hc_art_init(hd, cat, a, s) = intermediate.children.hc_art_init(hd, cat, a, s);
              }
            }// end ss.hTS
          }
        }// end hcTT
      }// end ss.hc1DS
    }// end a
  }// end ss.NS


  for(int s = 0; s <ss.NS; ++s) {
    for(int a = 0; a < pars.base.options.p_idx_fertility_first; ++a) {
      for(int hd = 0; hd < hc_ss.hc1DS; ++hd) {
        for(int cat = 0; cat < hc_ss.hcTT; ++cat) {
          intermediate.children.hc_art_init_total += intermediate.children.hc_art_init(hd, cat, a, s);
        }// end hcTT
      }// end ss.hc1DS
    }// end a
  }// end ss.NS

  //!!! TODO: fix order of for loop
  for(int s = 0; s <ss.NS; ++s) {
    for(int hd = 0; hd < hc_ss.hc1DS; ++hd) {
      for(int a = 0; a < pars.base.options.p_idx_fertility_first; ++a) {
        intermediate.children.hc_death_rate = 0.0;
        intermediate.children.hc_art_grad(0, hd, a, s) = 0.0;

        if(a < hc_ss.hc2_agestart){
          if(state_next.children.hc1_art_pop(0, hd, a, s) >0){
            intermediate.children.hc_death_rate =  cpars.hc_art_mort_rr(0, a, time_step) * 0.5 * (cpars.hc1_art_mort(hd, 0, a) + cpars.hc1_art_mort(hd, 1, a)) ;
            state_next.children.hc1_art_aids_deaths(0,hd, a, s) =  intermediate.children.hc_death_rate * state_next.children.hc1_art_pop(0, hd, a, s)  ;
            intermediate.children.hc_art_grad(0,hd, a, s) -= state_next.children.hc1_art_aids_deaths(0,hd, a, s) ;
            state_next.children.hc1_art_pop(0, hd,  a, s) += intermediate.children.hc_art_grad(0, hd, a, s) ;
          }
        }else{
          if(hd < hc_ss.hc2DS){
            if(state_next.children.hc2_art_pop(0, hd, a-hc_ss.hc2_agestart, s) >0){
            intermediate.children.hc_death_rate =  cpars.hc_art_mort_rr(0, a, time_step) * 0.5 * (cpars.hc2_art_mort(hd, 0, a-hc_ss.hc2_agestart) + cpars.hc2_art_mort(hd, 1, a-hc_ss.hc2_agestart));
            state_next.children.hc2_art_aids_deaths(0,hd, a-hc_ss.hc2_agestart, s) =  intermediate.children.hc_death_rate * state_next.children.hc2_art_pop(0, hd, a-hc_ss.hc2_agestart, s)  ;
            intermediate.children.hc_art_grad(0, hd, a, s) -= state_next.children.hc2_art_aids_deaths(0,hd, a-hc_ss.hc2_agestart, s) ;
            state_next.children.hc2_art_pop(0, hd,  a-hc_ss.hc2_agestart, s) += intermediate.children.hc_art_grad(0, hd, a, s) ;
          }
         }
        }
      }// end a
    }// end ss.hc1DS
  }// end ss.NS


  for(int s = 0; s <ss.NS; ++s) {
    for(int hd = 0; hd < hc_ss.hc1DS; ++hd) {
      for(int a = 0; a < pars.base.options.p_idx_fertility_first; ++a) {
        intermediate.children.hc_death_rate = 0;
        intermediate.children.hc_art_grad(2, hd, a, s) = 0.0;

        if(a < hc_ss.hc2_agestart){
          intermediate.children.hc_death_rate = cpars.hc_art_mort_rr(2, a, time_step) * cpars.hc1_art_mort(hd, 2, a);
          state_next.children.hc1_art_aids_deaths(2,hd, a, s) =  state_next.children.hc1_art_pop(2, hd, a, s) * intermediate.children.hc_death_rate;
          intermediate.children.hc_art_grad(2,hd, a, s) -= state_next.children.hc1_art_aids_deaths(2,hd, a, s) ;
          state_next.children.hc1_art_pop(2, hd,  a, s) += intermediate.children.hc_art_grad(2, hd, a, s) ;
        }else{
          if(hd < (hc_ss.hc2DS)){
            intermediate.children.hc_death_rate = cpars.hc_art_mort_rr(2, a, time_step) * cpars.hc2_art_mort(hd, 2, a-hc_ss.hc2_agestart);
            state_next.children.hc2_art_aids_deaths(2,hd, a-hc_ss.hc2_agestart, s) =  state_next.children.hc2_art_pop(2, hd, a-hc_ss.hc2_agestart, s) * intermediate.children.hc_death_rate;
            intermediate.children.hc_art_grad(2,hd, a, s) -= state_next.children.hc2_art_aids_deaths(2,hd, a-hc_ss.hc2_agestart, s) ;
            state_next.children.hc2_art_pop(2, hd,  a-hc_ss.hc2_agestart, s) += intermediate.children.hc_art_grad(2, hd, a, s) ;
          }
        }
      }// end a
    }// end ss.hc1DS
  }// end ss.NS

  //Progress ART to the correct time on ART
  //!!! TODO: fix order of for loop
  for(int hd = 0; hd < hc_ss.hc1DS; ++hd) {
    for(int a = 0; a < pars.base.options.p_idx_fertility_first; ++a) {
      for(int s = 0; s <ss.NS; ++s) {
        if(a < hc_ss.hc2_agestart){
          if(state_next.children.hc1_art_pop(0, hd, a, s) > 0){
            state_next.children.hc1_art_pop(1, hd, a, s) += state_next.children.hc1_art_pop(0, hd, a, s);
            state_next.children.hc1_art_pop(0, hd, a, s) -= state_next.children.hc1_art_pop(0, hd, a, s);
          }
        }else{
        //  if(state_next.children.hc2_art_pop(0, hd, a-ss.hc2_agestart, s) > 0){
            if(hd < (hc_ss.hc2DS)){
              state_next.children.hc2_art_pop(1, hd, a-hc_ss.hc2_agestart, s) += state_next.children.hc2_art_pop(0, hd, a-hc_ss.hc2_agestart, s);
              state_next.children.hc2_art_pop(0, hd, a-hc_ss.hc2_agestart, s) -= state_next.children.hc2_art_pop(0, hd, a-hc_ss.hc2_agestart, s);
            }
        //  }
        }
      }//end ss.NS
    }// end a
  }// end ss.hc1DS


  if ( !cpars.hc_art_isperc(time_step) && !cpars.hc_art_isperc(time_step-1) ){ // both numbers
    //Remove how many that are already on ART
    state_next.children.hc_art_num =  (cpars.hc_art_val(time_step) + cpars.hc_art_val(time_step-1)) / 2 ;
    for(int s = 0; s <ss.NS; ++s) {
      for(int a = 0; a < pars.base.options.p_idx_fertility_first; ++a) {
        for(int hd = 0; hd < hc_ss.hc1DS; ++hd) {
          for(int dur = 0; dur < ss.hTS; ++dur) {
            if(a < hc_ss.hc2_agestart){
              state_next.children.hc_art_num -= state_next.children.hc1_art_pop(dur, hd, a, s)  ;
            }else{
              if(hd < (hc_ss.hc2DS)){
                state_next.children.hc_art_num -= state_next.children.hc2_art_pop(dur, hd, a-hc_ss.hc2_agestart, s)   ;
              }
            }
          }// end ss.hTS
        }// end ss.hc1DS
      }// end a
    }// end ss.NS

    if(intermediate.children.hc_art_init_total < state_next.children.hc_art_num){
      state_next.children.hc_art_num = intermediate.children.hc_art_init_total;
    }
    if(state_next.children.hc_art_num < 0){
      state_next.children.hc_art_num =  0;
    }

  } else if (cpars.hc_art_isperc(time_step) && cpars.hc_art_isperc(time_step-1)){ // both percentages
    for(int s = 0; s <ss.NS; ++s) {
      for(int a = 0; a < pars.base.options.p_idx_fertility_first; ++a) {
        for(int hd = 0; hd < hc_ss.hc1DS; ++hd) {
          for(int dur = 0; dur < ss.hTS; ++dur) {
            if(a < hc_ss.hc2_agestart){
              state_next.children.hc_art_num += state_next.children.hc1_art_pop(dur, hd, a, s)  ;
              state_next.children.hc_art_num += state_next.children.hc1_art_aids_deaths(dur,hd, a, s) ;
            }else{
              if(hd < (hc_ss.hc2DS )){
                state_next.children.hc_art_num += state_next.children.hc2_art_pop(dur, hd, a-hc_ss.hc2_agestart, s)  ;
                state_next.children.hc_art_num += state_next.children.hc2_art_aids_deaths(dur,hd, a-hc_ss.hc2_agestart, s) ;
              }
            }
          } //end ss.hTS
        } // end ss.hC1_disease_stages
      } // end a
    } // end ss.NS

    state_next.children.hc_art_num =  state_next.children.hc_art_num * (cpars.hc_art_val(time_step) + cpars.hc_art_val(time_step-1)) / 2 ;

    //Remove how many that are already on ART
    for(int s = 0; s <ss.NS; ++s) {
      for(int a = 0; a < pars.base.options.p_idx_fertility_first; ++a) {
        for(int hd = 0; hd < hc_ss.hc1DS; ++hd) {
          for(int dur = 0; dur < ss.hTS; ++dur) {
            if(a < hc_ss.hc2_agestart){
              state_next.children.hc_art_num -= state_next.children.hc1_art_pop(dur, hd, a, s)  ;
            }else{
              if(hd < (hc_ss.hc2DS)){
                state_next.children.hc_art_num -= state_next.children.hc2_art_pop(dur, hd, a-hc_ss.hc2_agestart, s)  ;
              }
            }
          } // end ss.hTS
        } // end ss.hC1_disease_stages
      } // end a
    } // end ss.NS
    if(intermediate.children.hc_art_init_total < state_next.children.hc_art_num){
      state_next.children.hc_art_num = intermediate.children.hc_art_init_total;
    }



  } else if (cpars.hc_art_isperc(time_step) && !cpars.hc_art_isperc(time_step-1)){ // num to percentage
    //Remove how many that are already on ART
    for(int s = 0; s <ss.NS; ++s) {
      for(int a = 0; a < pars.base.options.p_idx_fertility_first; ++a) {
        for(int hd = 0; hd < hc_ss.hc1DS; ++hd) {
          for(int dur = 0; dur < ss.hTS; ++dur) {
            if(a < hc_ss.hc2_agestart){
              state_next.children.hc_art_num += state_next.children.hc1_art_pop(dur, hd, a, s) +  state_next.children.hc1_art_aids_deaths(dur,hd, a, s)  ;
            }else{
              if(hd < (hc_ss.hc2DS)){
                state_next.children.hc_art_num += state_next.children.hc2_art_pop(dur, hd, a-hc_ss.hc2_agestart, s) +  state_next.children.hc2_art_aids_deaths(dur,hd, a-hc_ss.hc2_agestart, s)  ;
              }
            }
          } // end ss.hTS
        } // end ss.hc1DS
      } //end a
    } //end ss.NS
    state_next.children.hc_art_num = (cpars.hc_art_val(time_step-1) + (state_next.children.hc_art_num * cpars.hc_art_val(time_step))) / 2 ;

    for(int s = 0; s <ss.NS; ++s) {
      for(int a = 0; a < pars.base.options.p_idx_fertility_first; ++a) {
        for(int hd = 0; hd < hc_ss.hc1DS; ++hd) {
          for(int dur = 0; dur < ss.hTS; ++dur) {
            if(a < hc_ss.hc2_agestart){
              state_next.children.hc_art_num -= state_next.children.hc1_art_pop(dur, hd, a, s);
            }else{
              if(hd < (hc_ss.hc2DS )){
                state_next.children.hc_art_num -= state_next.children.hc2_art_pop(dur, hd, a-hc_ss.hc2_agestart, s);
              }
            }
          } // end ss.hTS
        } // end ss.hc1DS
      } // end a
    } // end ss.NS

    if(state_next.children.hc_art_num < 0){
      state_next.children.hc_art_num = 0;
    }
    if(intermediate.children.hc_art_init_total < state_next.children.hc_art_num){
      state_next.children.hc_art_num = intermediate.children.hc_art_init_total;
    }

  } else if (cpars.hc_art_isperc(time_step-1) && !cpars.hc_art_isperc(time_step)){ //percentage to num


    for(int s = 0; s <ss.NS; ++s) {
      for(int a = 0; a < pars.base.options.p_idx_fertility_first; ++a) {
        for(int hd = 0; hd < hc_ss.hc1DS; ++hd) {
          for(int dur = 0; dur < ss.hTS; ++dur) {
            if(a < hc_ss.hc2_agestart){
              state_next.children.hc_art_num -= state_next.children.hc1_art_pop(dur, hd, a, s)  ;
            }else{
              if(hd < (hc_ss.hc2DS)){
                state_next.children.hc_art_num -= state_next.children.hc2_art_pop(dur, hd, a-hc_ss.hc2_agestart, s)   ;
              }
            }
          } // end ss.hTS
        } //end ss.hc1DS
      } // end a
    } //end ss.NS


   state_next.children.hc_art_num = (state_curr.children.hc_art_num + cpars.hc_art_val(time_step)) / 2 ;

    //Remove how many that are already on ART
    for(int s = 0; s <ss.NS; ++s) {
      for(int a = 0; a < pars.base.options.p_idx_fertility_first; ++a) {
        for(int hd = 0; hd < hc_ss.hc1DS; ++hd) {
          for(int dur = 0; dur < ss.hTS; ++dur) {
            if(a < hc_ss.hc2_agestart){
              state_next.children.hc_art_num -= state_next.children.hc1_art_pop(dur, hd, a, s);
            }else{
              if(hd < (hc_ss.hc2DS)){
                state_next.children.hc_art_num -= state_next.children.hc2_art_pop(dur, hd, a-hc_ss.hc2_agestart, s);
              }
            }
          } //end ss.hTS
        } //end ss.hc1DS
      } //end a
    } //end ss.NS

    if(state_next.children.hc_art_num < 0){
      state_next.children.hc_art_num = 0;
    }
    if(intermediate.children.hc_art_init_total < state_next.children.hc_art_num){
      state_next.children.hc_art_num = intermediate.children.hc_art_init_total;
    }
  }



  for(int s = 0; s <ss.NS; ++s) {
    for(int a = 0; a < pars.base.options.p_idx_fertility_first; ++a) {
      for(int hd = 0; hd < hc_ss.hc1DS; ++hd) {
        for(int cat = 0; cat < hc_ss.hcTT; ++cat) {
          intermediate.children.hc_initByAge +=  intermediate.children.hc_art_init(hd, cat, a, s) * cpars.hc_art_init_dist(a, time_step);
        } //end hcTT
      } // end ss.hc1DS
    } //end a
  } // end ss.NS

  if(intermediate.children.hc_initByAge == 0.0){
    intermediate.children.hc_adj = 1.0 ;
  }else{
    intermediate.children.hc_adj = state_next.children.hc_art_num / intermediate.children.hc_initByAge ;
  }
  for(int s = 0; s <ss.NS; ++s) {
    for(int cat = 0; cat < hc_ss.hcTT; ++cat) {
      for(int a = 0; a < pars.base.options.p_idx_fertility_first; ++a) {
        for(int hd = 0; hd < hc_ss.hc1DS; ++hd) {
          if((intermediate.children.hc_adj * cpars.hc_art_init_dist(a, time_step)) > 1.0){
            intermediate.children.hc_art_scalar = 1.0;
          }else{
            intermediate.children.hc_art_scalar = intermediate.children.hc_adj * cpars.hc_art_init_dist(a, time_step);
          }
          if(state_next.children.hc_art_num > 0.0){
            intermediate.children.hc_art_scalar = intermediate.children.hc_art_scalar;
          }else{
            intermediate.children.hc_art_scalar =  0.0;
          }
          if(a < hc_ss.hc2_agestart){
            state_next.children.hc1_art_pop(0, hd, a, s) +=  intermediate.children.hc_art_scalar * intermediate.children.hc_art_init(hd, cat, a, s) ;
          }else{
            if(hd < (hc_ss.hc2DS)){
              state_next.children.hc2_art_pop(0, hd, a - hc_ss.hc2_agestart, s) +=  intermediate.children.hc_art_scalar * intermediate.children.hc_art_init(hd, cat, a, s) ;
            }
          }
          if(a < hc_ss.hc2_agestart){
            state_next.children.hc1_hiv_pop(hd, cat, a, s) -=  intermediate.children.hc_art_scalar * intermediate.children.hc_art_init(hd, cat, a, s) ;
          }else{
            if(hd < (hc_ss.hc2DS )){
              state_next.children.hc2_hiv_pop(hd, cat, a - hc_ss.hc2_agestart, s) -=  intermediate.children.hc_art_scalar * intermediate.children.hc_art_init(hd, cat, a, s) ;
            }
          }

        } //end ss.hc1DS
      } // end a
    } // end hcTT
  } // end ss.NS
}


template<typename ModelVariant, typename real_type>
void run_child_art_mortality(int time_step,
                             const Parameters<ModelVariant, real_type> &pars,
                             const State<ModelVariant, real_type> &state_curr,
                             State<ModelVariant, real_type> &state_next,
                             IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "run_hiv_child_infections can only be called for model variants where run_child_model is true");
  constexpr auto ss = StateSpace<ModelVariant>().base;
  constexpr auto hc_ss = StateSpace<ModelVariant>().children;
  const auto cpars = pars.children.children;


  for(int hd = 0; hd < hc_ss.hc1DS; ++hd) {
    for(int a = 0; a < pars.base.options.p_idx_fertility_first; ++a) {
      for(int s = 0; s <ss.NS; ++s) {
        intermediate.children.hc_death_rate = 0.0;
        intermediate.children.hc_art_grad(0, hd, a, s) = 0.0;
        if(a < hc_ss.hc2_agestart){
          if(state_next.children.hc1_art_pop(0, hd, a, s) > 0){
            intermediate.children.hc_death_rate = cpars.hc_art_mort_rr(0, a, time_step) * 0.5 * (cpars.hc1_art_mort(hd, 0, a) + cpars.hc1_art_mort(hd, 1, a));
            state_next.children.hc1_art_aids_deaths(0, hd, a, s) = intermediate.children.hc_death_rate * state_next.children.hc1_art_pop(0, hd, a, s);
            intermediate.children.hc_art_grad(0, hd, a, s) -= state_next.children.hc1_art_aids_deaths(0, hd, a, s);
            state_next.children.hc1_art_pop(0, hd, a, s) += intermediate.children.hc_art_grad(0, hd, a, s);
          }
        }else{
          if(hd < hc_ss.hc2DS){
            intermediate.children.hc_death_rate = cpars.hc_art_mort_rr(0, a, time_step) * 0.5 * (cpars.hc2_art_mort(hd, 0, a-hc_ss.hc2_agestart) + cpars.hc2_art_mort(hd, 1, a-hc_ss.hc2_agestart));
            state_next.children.hc2_art_aids_deaths(0, hd, a-hc_ss.hc2_agestart, s)  = intermediate.children.hc_death_rate * state_next.children.hc2_art_pop(0, hd, a-hc_ss.hc2_agestart, s);
            intermediate.children.hc_art_grad(0, hd, a, s) -= state_next.children.hc2_art_aids_deaths(0, hd, a-hc_ss.hc2_agestart, s);
            state_next.children.hc2_art_pop(0, hd, a-hc_ss.hc2_agestart, s) += intermediate.children.hc_art_grad(0, hd, a, s);
          }
        }
      }//end ss.NS
    }// end a
  }// end ss.hc1DS

  for(int hd = 0; hd < hc_ss.hc1DS; ++hd) {
    for(int a = 0; a < pars.base.options.p_idx_fertility_first; ++a) {
      for(int s = 0; s <ss.NS; ++s) {
        if(a < hc_ss.hc2_agestart){
          if(state_next.children.hc1_art_pop(1, hd, a, s) > 0){
            state_next.children.hc1_art_pop(2, hd, a, s) += state_next.children.hc1_art_pop(1, hd, a, s);
            state_next.children.hc1_art_pop(1, hd, a, s) -= state_next.children.hc1_art_pop(1, hd, a, s);
          }
        }else{
          if(state_next.children.hc2_art_pop(1, hd, a-hc_ss.hc2_agestart, s) > 0){
            if(hd < (hc_ss.hc2DS)){
              state_next.children.hc2_art_pop(2, hd, a-hc_ss.hc2_agestart, s) += state_next.children.hc2_art_pop(1, hd, a-hc_ss.hc2_agestart, s);
              state_next.children.hc2_art_pop(1, hd, a-hc_ss.hc2_agestart, s) -= state_next.children.hc2_art_pop(1, hd, a-hc_ss.hc2_agestart, s);
            }
          }
        }
      }//end ss.NS
    }// end a
  }// end ss.hc1DS
}


}// namespace internal

template<typename ModelVariant, typename real_type>
void run_child_model_simulation(int time_step,
                                const Parameters<ModelVariant, real_type> &pars,
                                const State<ModelVariant, real_type> &state_curr,
                                State<ModelVariant, real_type> &state_next,
                                internal::IntermediateData<ModelVariant, real_type> &intermediate) {
 internal::run_child_hiv_infections(time_step, pars, state_curr, state_next, intermediate);
 internal::run_child_natural_history(time_step, pars, state_curr, state_next, intermediate);
 internal::run_child_hiv_mort(time_step, pars, state_curr, state_next, intermediate);
 internal::add_child_grad(time_step, pars, state_curr, state_next, intermediate);
 // this function may need to be broken up, its around 350 lines
 // !!!TODO: also need to fix the looping order for some loops
 // !!!TODO: put this in an if statement to only run if the first year of ART has passed
 internal::run_child_art_initiation(time_step, pars, state_curr, state_next, intermediate);
 internal::run_child_art_mortality(time_step, pars, state_curr, state_next, intermediate);

}

} // namespace leapfrog
