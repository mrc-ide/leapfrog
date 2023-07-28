#pragma once

#include "types.hpp"

namespace leapfrog {

namespace internal {

template<HivAgeStratification S, typename real_type>
void run_child_ageing(int time_step,
                     const Parameters<real_type> &pars,
                     const State<S, real_type> &state_curr,
                     State<S, real_type> &state_next,
                     IntermediateData<S, real_type> &intermediate) {
  const auto demog = pars.demography;
  const auto cpars = pars.children;

  constexpr auto ss = StateSpace<S>();
  for (int g = 0; g < ss.num_genders; ++g) {
    //less than 5 because there is a cd4 transition between ages 4 and 5
    for (int af = 1; af < ss.hc2_agestart; ++af) {
      for (int hm = 0; hm < ss.hc1_disease_stages; ++hm) {
        for (int cat = 0 ; cat < ss.hTM; ++cat) {
          state_next.hc1_hiv_pop(hm, cat, af, g) += state_curr.hc1_hiv_pop(hm, cat, af-1, g) * demog.survival(af, g, time_step);

        }
        for (int dur = 0; dur < ss.treatment_stages; ++dur) {
          state_next.hc1_art_pop(dur, hm, af, g) += state_curr.hc1_art_pop(dur, hm, af-1, g) * demog.survival(af, g, time_step);
        }

      }
    }
  }

  for (int g = 0; g < ss.num_genders; ++g) {
    for (int hm = 0; hm < ss.hc1_disease_stages; ++hm) {
      for (int hm_alt = 0; hm_alt < ss.hc2_disease_stages; ++hm_alt) {
        for (int cat = 0 ; cat < ss.hTM; ++cat) {
          state_next.hc2_hiv_pop(hm_alt, cat, 0, g) +=  state_curr.hc1_hiv_pop(hm, cat, ss.hc1_ageend, g) * demog.survival(ss.hc2_agestart, g, time_step) * cpars.hc_cd4_transition(hm_alt, hm);
        }
        for (int dur = 0; dur < ss.treatment_stages; ++dur) {
          state_next.hc2_art_pop(dur, hm_alt, 0, g) += state_curr.hc1_art_pop(dur, hm, ss.hc1_ageend, g) * demog.survival(ss.hc2_agestart, g, time_step) * cpars.hc_cd4_transition(hm_alt, hm);
        }
      }
    }
  }

  for (int g = 0; g < ss.num_genders; ++g) {
    for (int af = (ss.hc2_agestart + 1); af < pars.options.fertility_first_age_group; ++af) {
      for (int hm = 0; hm < ss.hc2_disease_stages; ++hm) {
        for (int cat = 0 ; cat < ss.hTM; ++cat) {
          state_next.hc2_hiv_pop(hm, cat, af - ss.hc2_agestart, g) += state_curr.hc2_hiv_pop(hm, cat, af- ss.hc2_agestart-1, g) * demog.survival(af, g, time_step);

        }
        for (int dur = 0; dur < ss.treatment_stages; ++dur) {
          state_next.hc2_art_pop(dur, hm, af - ss.hc2_agestart, g) += state_curr.hc2_art_pop(dur, hm, af- ss.hc2_agestart-1, g) * demog.survival(af, g, time_step);
        }

      }
    }
  }



}

template<HivAgeStratification S, typename real_type>
void run_child_hiv_infections(int time_step,
                              const Parameters<real_type> &pars,
                              const State<S, real_type> &state_curr,
                              State<S, real_type> &state_next,
                              IntermediateData<S, real_type> &intermediate) {
  constexpr auto ss = StateSpace<S>();
  const auto cpars = pars.children;

  for (int g = 0; g < ss.num_genders; ++g) {
    // Run only first 5 age groups in total population 0, 1, 2, 3, 4
    for (int af = 0; af < ss.hc2_agestart; ++af) {
      if (cpars.hc_nosocomial(time_step) > 0) {
        // Divide by 10 because we want to evenly distribute over 2 genders and 5 age groups
        state_next.infections(af, g) = cpars.hc_nosocomial(time_step) / (5.0 * ss.num_genders);
        state_next.hiv_population(af, g) += state_next.infections(af, g);

        for (int hm = 0; hm < ss.hc1_disease_stages; ++hm) {
          // putting them all in perinatal hTM to match spec nosocomial
          if (cpars.hc1_cd4_dist(hm) > 0) {
            state_next.hc1_hiv_pop(hm, 0, af, g) += state_next.infections(af, g) * cpars.hc1_cd4_dist(hm);
          }
        }
      }
    }
  }
}


template<HivAgeStratification S, typename real_type>
void run_child_natural_history(int time_step,
                              const Parameters<real_type> &pars,
                              const State<S, real_type> &state_curr,
                              State<S, real_type> &state_next,
                              IntermediateData<S, real_type> &intermediate) {
  constexpr auto ss = StateSpace<S>();
  const auto cpars = pars.children;


 for (int g = 0; g < ss.num_genders; ++g) {
   for (int hm = 0; hm < ss.hc1_disease_stages; ++hm) {
     for (int af = 0; af < ss.hc2_agestart; ++af) {
       for (int cat = 0; cat < ss.hTM; ++cat) {
         intermediate.hc_posthivmort(hm, cat, af, g) += state_next.hc1_hiv_pop(hm, cat, af, g) - (1.0 - cpars.ctx_effect * cpars.ctx_val(time_step)) * state_curr.hc1_hiv_pop(hm, cat, af, g) * cpars.hc1_cd4_mort(hm, cat, af);
         state_next.hc1_noart_aids_deaths(hm, cat, af, g) += (1 - cpars.ctx_effect * cpars.ctx_val(time_step)) * state_curr.hc1_hiv_pop(hm, cat, af, g) * cpars.hc1_cd4_mort(hm, cat, af)  ;
       }
     }
   }
 }


 for (int g = 0; g < ss.num_genders; ++g) {
   for (int hm = 0; hm < ss.hc2_disease_stages; ++hm) {
     for (int af = ss.hc2_agestart; af < pars.options.fertility_first_age_group; ++af) {
       for (int cat = 0; cat < ss.hTM; ++cat) {
         intermediate.hc_posthivmort(hm, cat, af, g) += state_next.hc2_hiv_pop(hm, cat, af - ss.hc2_agestart, g) - (1 - cpars.ctx_effect * cpars.ctx_val(time_step)) * state_curr.hc2_hiv_pop(hm, cat, af - ss.hc2_agestart, g) * cpars.hc2_cd4_mort(hm, cat, af - ss.hc2_agestart);
         state_next.hc2_noart_aids_deaths(hm, cat, af, g) +=  (1 - cpars.ctx_effect * cpars.ctx_val(time_step)) * state_curr.hc2_hiv_pop(hm, cat, af - ss.hc2_agestart, g) * cpars.hc2_cd4_mort(hm, cat, af - ss.hc2_agestart); // output hiv deaths, aggregated across transmission category
       }
     }
   }
 }

 //progress through CD4 categories
 for (int g = 0; g < ss.num_genders; ++g) {
   for (int hm = 1; hm < ss.hc1_disease_stages; ++hm) {
     for (int af = 0; af < ss.hc2_agestart; ++af) {
       for (int cat = 0; cat < ss.hTM; ++cat) {
         intermediate.hc_grad(hm - 1, cat, af, g) -=  (intermediate.hc_posthivmort(hm - 1, cat, af, g) * cpars.hc1_cd4_prog(hm - 1) + state_next.hc1_hiv_pop(hm - 1, cat, af, g) * cpars.hc1_cd4_prog(hm - 1)) / 2; //moving to next cd4 category
         intermediate.hc_grad(hm, cat, af, g) += (intermediate.hc_posthivmort(hm - 1, cat, af, g) * cpars.hc1_cd4_prog(hm - 1) + state_next.hc1_hiv_pop(hm - 1, cat, af, g) * cpars.hc1_cd4_prog(hm - 1)) / 2; //moving into this cd4 category
       }
     }
   }
 }

 //progress through CD4 categories
 for (int g = 0; g < ss.num_genders; ++g) {
   for (int hm = 1; hm < ss.hc2_disease_stages; ++hm) {
     for (int af = ss.hc2_agestart; af < pars.options.fertility_first_age_group; ++af) {
       for (int cat = 0; cat < ss.hTM; ++cat) {
         intermediate.hc_grad(hm - 1, cat, af, g) -= (intermediate.hc_posthivmort(hm - 1, cat, af, g) * cpars.hc2_cd4_prog(hm - 1) + state_next.hc2_hiv_pop(hm - 1, cat, af - ss.hc2_agestart, g) * cpars.hc2_cd4_prog(hm - 1)) / 2; //moving to next cd4 category
         intermediate.hc_grad(hm, cat, af, g) += (intermediate.hc_posthivmort(hm - 1, cat, af, g) * cpars.hc2_cd4_prog(hm - 1) + state_next.hc2_hiv_pop(hm - 1, cat, af - ss.hc2_agestart, g) * cpars.hc2_cd4_prog(hm - 1)) / 2; //moving into this cd4 category
       }
     }
   }
 }
}

template<HivAgeStratification S, typename real_type>
void run_child_hiv_mort(int time_step,
                               const Parameters<real_type> &pars,
                               const State<S, real_type> &state_curr,
                               State<S, real_type> &state_next,
                               IntermediateData<S, real_type> &intermediate) {
  constexpr auto ss = StateSpace<S>();
  const auto cpars = pars.children;


  for (int g = 0; g < ss.num_genders; ++g) {
    for (int hm = 0; hm < ss.hc2_disease_stages; ++hm) {
      for (int af = 0; af < pars.options.fertility_first_age_group; ++af) {
        for (int cat = 0; cat < ss.hTM; ++cat) {
          if(af < ss.hc2_agestart){
            intermediate.hc_grad(hm, cat, af, g) -= (1 - cpars.ctx_effect * cpars.ctx_val(time_step)) * state_next.hc1_hiv_pop(hm, cat, af, g) * cpars.hc1_cd4_mort(hm, cat, af)  ;
          }else{
            intermediate.hc_grad(hm, cat, af, g) -= (1 - cpars.ctx_effect * cpars.ctx_val(time_step)) * state_next.hc2_hiv_pop(hm, cat, af - ss.hc2_agestart, g) * cpars.hc2_cd4_mort(hm, cat, af - ss.hc2_agestart)  ;

          }
        }
      }
    }
  }
}

template<HivAgeStratification S, typename real_type>
void add_child_grad(int time_step,
                         const Parameters<real_type> &pars,
                         const State<S, real_type> &state_curr,
                         State<S, real_type> &state_next,
                         IntermediateData<S, real_type> &intermediate) {
  constexpr auto ss = StateSpace<S>();
  const auto cpars = pars.children;

  //add on transitions
  for (int g = 0; g < ss.num_genders; ++g) {
    for (int hm = 0; hm < ss.hc1_disease_stages; ++hm) {
      for (int af = 0; af < pars.options.fertility_first_age_group; ++af) {
        for (int cat = 0; cat < ss.hTM; ++cat) {
          if(af < ss.hc2_agestart){
            state_next.hc1_hiv_pop(hm, cat, af, g) +=intermediate.hc_grad(hm, cat, af, g);
          }else{
          state_next.hc2_hiv_pop(hm, cat, af- ss.hc2_agestart, g) +=intermediate.hc_grad(hm, cat, af, g);
        }
      }
    }
  }

}

}

template<HivAgeStratification S, typename real_type>
void run_child_art_initiation(int time_step,
                    const Parameters<real_type> &pars,
                    const State<S, real_type> &state_curr,
                    State<S, real_type> &state_next,
                    IntermediateData<S, real_type> &intermediate) {
  constexpr auto ss = StateSpace<S>();
  const auto cpars = pars.children;

  //all children under a certain age eligible for ART
  for(int g = 0; g < ss.num_genders; ++g)  {
    for(int cat = 0; cat < ss.hTM; ++cat) {
      for(int af = 0; af < pars.options.fertility_first_age_group; ++af) {
        for(int hm = 0; hm < ss.hc1_disease_stages; ++hm) {

          if(af >  (cpars.hc_art_elig_cd4(af, time_step) - 2)){
            if(af < ss.hc2_agestart){
              intermediate.hc_art_need(hm, cat, af, g) += state_next.hc1_hiv_pop(hm, cat, af, g);
            }else{
              intermediate.hc_art_need(hm, cat, af, g) += state_next.hc2_hiv_pop(hm, cat, af - ss.hc2_agestart, g);
            }
            }

            if(af >  cpars.hc_art_elig_age(time_step) & hm > (cpars.hc_art_elig_cd4(af, time_step) - 2)){
              if(af < ss.hc2_agestart){
                intermediate.hc_art_need(hm, cat, af, g) += state_next.hc1_hiv_pop(hm, cat, af, g);
              }else{
                intermediate.hc_art_need(hm, cat, af, g) += state_next.hc2_hiv_pop(hm, cat, af, g);
              }
            }

        } // end ss.hc1_disease_stages
      } // end af
    } // end hTM
  } // end ss.num_genders

  for(int g = 0; g < ss.num_genders; ++g) {
    for(int af = 0; af < pars.options.fertility_first_age_group; ++af) {
      for(int hm = 0; hm < ss.hc1_disease_stages; ++hm) {
        for(int cat = 0; cat < ss.hTM; ++cat) {
          intermediate.hc_art_num += intermediate.hc_art_need(hm, cat, af, g);
        } // end hTM
      } // end ss.hc1_disease_stages
    } // end af
  } // end ss.num_genders

  //how many should initialize ART
  for(int g = 0; g < ss.num_genders; ++g) {
    for(int af = 0; af < pars.options.fertility_first_age_group; ++af) {
      for(int hm = 0; hm < ss.hc1_disease_stages; ++hm) {
        for(int cat = 0; cat < ss.hTM; ++cat) {
          if(cpars.hc_art_val(time_step) > 0){
            intermediate.hc_art_init(hm, cat, af, g) += intermediate.hc_art_need(hm, cat, af, g);
            for(int dur = 0; dur < ss.treatment_stages; ++dur) {
              if(intermediate.hc_art_init(hm, cat, af, g) < 0){
                intermediate.hc_art_init(hm, cat, af, g) = 0.0;
              }else{
                intermediate.hc_art_init(hm, cat, af, g) = intermediate.hc_art_init(hm, cat, af, g);
              }
            }// end ss.treatment_stages
          }
        }// end hTM
      }// end ss.hc1_disease_stages
    }// end af
  }// end ss.num_genders


  for(int g = 0; g < ss.num_genders; ++g) {
    for(int af = 0; af < pars.options.fertility_first_age_group; ++af) {
      for(int hm = 0; hm < ss.hc1_disease_stages; ++hm) {
        for(int cat = 0; cat < ss.hTM; ++cat) {
          intermediate.hc_art_init_total += intermediate.hc_art_init(hm, cat, af, g);
        }// end hTM
      }// end ss.hc1_disease_stages
    }// end af
  }// end ss.num_genders

  //!!! TODO: fix order of for loop
  for(int g = 0; g < ss.num_genders; ++g) {
    for(int hm = 0; hm < ss.hc1_disease_stages; ++hm) {
      for(int af = 0; af < pars.options.fertility_first_age_group; ++af) {
        if(af < ss.hc2_agestart){
          intermediate.hc_death_rate =  cpars.hc_art_mort_rr(0, af, time_step) * 0.5 * (cpars.hc1_art_mort(hm, 0, af) + cpars.hc1_art_mort(hm, 1, af)) ;
          state_next.hc1_art_aids_deaths(0,hm, af, g) =  intermediate.hc_death_rate * state_next.hc1_art_pop(0, hm, af, g)  ;
          intermediate.hc_art_grad(0,hm, af, g) -= state_next.hc1_art_aids_deaths(0,hm, af, g) ;
          state_next.hc1_art_pop(0, hm,  af, g) += intermediate.hc_art_grad(0, hm, af, g) ;
        }else{
          if(hm < (ss.hc2_disease_stages - 1)){
            intermediate.hc_death_rate =  cpars.hc_art_mort_rr(0, af, time_step) * 0.5 * (cpars.hc2_art_mort(hm, 0, af-ss.hc2_agestart) + cpars.hc2_art_mort(hm, 1, af-ss.hc2_agestart));
            state_next.hc2_art_aids_deaths(0,hm, af-ss.hc2_agestart, g) =  intermediate.hc_death_rate * state_next.hc2_art_pop(0, hm, af-ss.hc2_agestart, g)  ;
            intermediate.hc_art_grad(0, hm, af, g) -= state_next.hc2_art_aids_deaths(0,hm, af-ss.hc2_agestart, g) ;
            state_next.hc2_art_pop(0, hm,  af-ss.hc2_agestart, g) += intermediate.hc_art_grad(0, hm, af, g) ;
          }
        }
        intermediate.hc_art_grad(0, hm, af, g) = 0.0;
      }// end af
    }// end ss.hc1_disease_stages
  }// end ss.num_genders


  for(int g = 0; g < ss.num_genders; ++g) {
    for(int hm = 0; hm < ss.hc1_disease_stages; ++hm) {
      for(int af = 0; af < pars.options.fertility_first_age_group; ++af) {
        intermediate.hc_death_rate = 0;
        if(af < ss.hc2_agestart){
          intermediate.hc_death_rate = cpars.hc_art_mort_rr(2, af, time_step) * cpars.hc1_art_mort(hm, 2, af);
          state_next.hc1_art_aids_deaths(2,hm, af, g) =  state_next.hc1_art_pop(2, hm, af, g) * intermediate.hc_death_rate;
          intermediate.hc_art_grad(2,hm, af, g) -= state_next.hc1_art_aids_deaths(2,hm, af, g) ;
          state_next.hc1_art_pop(2, hm,  af, g) += intermediate.hc_art_grad(2, hm, af, g) ;
        }else{
          if(hm < (ss.hc2_disease_stages - 1)){
            intermediate.hc_death_rate = cpars.hc_art_mort_rr(2, af, time_step) * cpars.hc2_art_mort(hm, 2, af-ss.hc2_agestart);
            state_next.hc2_art_aids_deaths(2,hm, af-ss.hc2_agestart, g) =  state_next.hc2_art_pop(2, hm, af-ss.hc2_agestart, g) * intermediate.hc_death_rate;
            intermediate.hc_art_grad(2,hm, af, g) -= state_next.hc2_art_aids_deaths(2,hm, af-ss.hc2_agestart, g) ;
            state_next.hc2_art_pop(2, hm,  af-ss.hc2_agestart, g) += intermediate.hc_art_grad(2, hm, af, g) ;
          }
        }
      }// end af
    }// end ss.hc1_disease_stages
  }// end ss.num_genders

  //Progress ART to the correct time on ART
  //!!! TODO: fix order of for loop
  for(int hm = 0; hm < ss.hc1_disease_stages; ++hm) {
    for(int af = 0; af < pars.options.fertility_first_age_group; ++af) {
      for(int g = 0; g < ss.num_genders; ++g) {
        if(af < ss.hc2_agestart){
          if(state_next.hc1_art_pop(0, hm, af, g) > 0){
            state_next.hc1_art_pop(1, hm, af, g) += state_next.hc1_art_pop(0, hm, af, g);
            state_next.hc1_art_pop(0, hm, af, g) -= state_next.hc1_art_pop(0, hm, af, g);
          }
        }else{
          if(state_next.hc2_art_pop(0, hm, af-ss.hc2_agestart, g) > 0){
            if(hm < (ss.hc2_disease_stages - 1)){
              state_next.hc2_art_pop(1, hm, af-ss.hc2_agestart, g) += state_next.hc2_art_pop(0, hm, af-ss.hc2_agestart, g);
              state_next.hc2_art_pop(0, hm, af-ss.hc2_agestart, g) -= state_next.hc2_art_pop(0, hm, af-ss.hc2_agestart, g);
            }
          }
        }
      }//end ss.num_genders
    }// end af
  }// end ss.hc1_disease_stages


  if ( !cpars.hc_art_isperc(time_step) && !cpars.hc_art_isperc(time_step-1) ){ // both numbers
    //Remove how many that are already on ART
    intermediate.hc_art_num =  (cpars.hc_art_val(time_step) + cpars.hc_art_val(time_step-1)) / 2 ;
    for(int g = 0; g < ss.num_genders; ++g) {
      for(int af = 0; af < pars.options.fertility_first_age_group; ++af) {
        for(int hm = 0; hm < ss.hc1_disease_stages; ++hm) {
          for(int dur = 0; dur < ss.treatment_stages; ++dur) {
            if(af < ss.hc2_agestart){
              intermediate.hc_art_num -= state_next.hc1_art_pop(dur, hm, af, g)  ;
            }else{
              if(hm < (ss.hc2_disease_stages - 1)){
                intermediate.hc_art_num -= state_next.hc2_art_pop(dur, hm, af-ss.hc2_agestart, g)   ;
              }
            }
          }// end ss.treatment_stages
        }// end ss.hc1_disease_stages
      }// end af
    }// end ss.num_genders

    if(intermediate.hc_art_init_total < intermediate.hc_art_num){
      intermediate.hc_art_num = intermediate.hc_art_init_total;
    }
    if(intermediate.hc_art_num < 0){
      intermediate.hc_art_num =  0;
    }

  } else if (cpars.hc_art_isperc(time_step) && cpars.hc_art_isperc(time_step-1)){ // both percentages
    for(int g = 0; g < ss.num_genders; ++g) {
      for(int af = 0; af < pars.options.fertility_first_age_group; ++af) {
        for(int hm = 0; hm < ss.hc1_disease_stages; ++hm) {
          for(int dur = 0; dur < ss.treatment_stages; ++dur) {
            if(af < ss.hc2_agestart){
              intermediate.hc_art_num += state_next.hc1_art_pop(dur, hm, af, g)  ;
              intermediate.hc_art_num += state_next.hc1_art_aids_deaths(dur,hm, af, g) ;

            }else{
              if(hm < (ss.hc2_disease_stages - 1)){
                intermediate.hc_art_num += state_next.hc2_art_pop(dur, hm, af-ss.hc2_agestart, g)  ;
                intermediate.hc_art_num += state_next.hc2_art_aids_deaths(dur,hm, af-ss.hc2_agestart, g) ;

              }
            }
          } //end ss.treatment_stages
        } // end ss.hC1_disease_stages
      } // end af
    } // end ss.num_genders

    intermediate.hc_art_num =  intermediate.hc_art_num * (cpars.hc_art_val(time_step) + cpars.hc_art_val(time_step-1)) / 2 ;

    //Remove how many that are already on ART
    for(int g = 0; g < ss.num_genders; g) {
      for(int af = 0; af < pars.options.fertility_first_age_group; af) {
        for(int hm = 0; hm < ss.hc1_disease_stages; hm) {
          for(int dur = 0; dur < ss.treatment_stages; dur) {
            if(af < ss.hc2_agestart){
              intermediate.hc_art_num -= state_next.hc1_art_pop(dur, hm, af, g)  ;
            }else{
              if(hm < (ss.hc2_disease_stages - 1)){
                intermediate.hc_art_num -= state_next.hc2_art_pop(dur, hm, af-ss.hc2_agestart, g)  ;
              }
            }
          } // end ss.treatment_stages
        } // end ss.hC1_disease_stages
      } // end af
    } // end ss.num_genders
    if(intermediate.hc_art_init_total < intermediate.hc_art_num){
      intermediate.hc_art_num = intermediate.hc_art_init_total;
    }



  } else if (cpars.hc_art_isperc(time_step) && !cpars.hc_art_isperc(time_step-1)){ // num to percentage
    //Remove how many that are already on ART
    double temp ;
    temp = 0.0;
    for(int g = 0; g < ss.num_genders; ++g) {
      for(int af = 0; af < pars.options.fertility_first_age_group; ++af) {
        for(int hm = 0; hm < ss.hc1_disease_stages; ++hm) {
          for(int dur = 0; dur < ss.treatment_stages; ++dur) {
            if(af < ss.hc2_agestart){
              intermediate.hc_art_num += state_next.hc1_art_pop(dur, hm, af, g) +  state_next.hc1_art_aids_deaths(dur,hm, af, g)  ;
            }else{
              if(hm < (ss.hc2_disease_stages - 1)){
                intermediate.hc_art_num += state_next.hc2_art_pop(dur, hm, af-ss.hc2_agestart, g) +  state_next.hc2_art_aids_deaths(dur,hm, af-ss.hc2_agestart, g)  ;
              }
            }
          } // end ss.treatment_stages
        } // end ss.hc1_disease_stages
      } //end af
    } //end ss.num_genders
    intermediate.hc_art_num = (cpars.hc_art_val(time_step-1) + (intermediate.hc_art_num * cpars.hc_art_val(time_step))) / 2 ;

    for(int g = 0; g < ss.num_genders; ++g) {
      for(int af = 0; af < pars.options.fertility_first_age_group; ++af) {
        for(int hm = 0; hm < ss.hc1_disease_stages; ++hm) {
          for(int dur = 0; dur < ss.treatment_stages; ++dur) {
            if(af < ss.hc2_agestart){
              intermediate.hc_art_num -= state_next.hc1_art_pop(dur, hm, af, g);
            }else{
              if(hm < (ss.hc2_disease_stages - 1)){
                intermediate.hc_art_num -= state_next.hc2_art_pop(dur, hm, af-ss.hc2_agestart, g);
              }
            }
          } // end ss.treatment_stages
        } // end ss.hc1_disease_stages
      } // end af
    } // end ss.num_genders

    if(intermediate.hc_art_num < 0){
      intermediate.hc_art_num = 0;
    }
    if(intermediate.hc_art_init_total < intermediate.hc_art_num){
      intermediate.hc_art_num = intermediate.hc_art_init_total;
    }

  } else if (cpars.hc_art_isperc(time_step-1) && !cpars.hc_art_isperc(time_step)){ //percentage to num


    for(int g = 0; g < ss.num_genders; ++g) {
      for(int af = 0; af < pars.options.fertility_first_age_group; ++af) {
        for(int hm = 0; hm < ss.hc1_disease_stages; ++hm) {
          for(int dur = 0; dur < ss.treatment_stages; ++dur) {
            if(af < ss.hc2_agestart){
              intermediate.hc_art_num -= state_next.hc1_art_pop(dur, hm, af, g)  ;
            }else{
              if(hm < (ss.hc2_disease_stages - 1)){
                intermediate.hc_art_num -= state_next.hc2_art_pop(dur, hm, af-ss.hc2_agestart, g)   ;
              }
            }
          } // end ss.treatment_stages
        } //end ss.hc1_disease_stages
      } // end af
    } //end ss.num_genders


    //TODO !!! Actually not sure how to do this, bc intermediates shouldn't have any time steps...
   // intermediate.hc_art_num = (intermediate.hc_art_num(time_step-1) + cpars.hc_art_val(time_step)) / 2 ;


    //Remove how many that are already on ART
    for(int g = 0; g < ss.num_genders; ++g) {
      for(int af = 0; af < pars.options.fertility_first_age_group; ++af) {
        for(int hm = 0; hm < ss.hc1_disease_stages; ++hm) {
          for(int dur = 0; dur < ss.treatment_stages; ++dur) {
            if(af < ss.hc2_agestart){
              intermediate.hc_art_num -= state_next.hc1_art_pop(dur, hm, af, g);
            }else{
              if(hm < (ss.hc2_disease_stages - 1)){
                intermediate.hc_art_num -= state_next.hc2_art_pop(dur, hm, af-ss.hc2_agestart, g);
              }
            }
          } //end ss.treatment_stages
        } //end ss.hc1_disease_stages
      } //end af
    } //end ss.num_genders

    if(intermediate.hc_art_num < 0){
      intermediate.hc_art_num = 0;
    }
    if(intermediate.hc_art_init_total < intermediate.hc_art_num){
      intermediate.hc_art_num = intermediate.hc_art_init_total;
    }
  }


  double initByAge = 0.0;
  for(int g = 0; g < ss.num_genders; ++g) {
    for(int af = 0; af < pars.options.fertility_first_age_group; ++af) {
      for(int hm = 0; hm < ss.hc1_disease_stages; ++hm) {
        for(int cat = 0; cat < ss.hTM; ++cat) {
          initByAge +=  intermediate.hc_art_init(hm, cat, af, g) * cpars.hc_art_init_dist(af, time_step);
        } //end hTM
      } // end ss.hc1_disease_stages
    } //end af
  } // end ss.num_genders

  double adj ;
  if(initByAge == 0){
    adj = 1.0 ;
  }else{
    adj = intermediate.hc_art_num / initByAge ;
  }
  for(int g = 0; g < ss.num_genders; ++g) {
    for(int cat = 0; cat < ss.hTM; ++cat) {
      for(int af = 0; af < pars.options.fertility_first_age_group; ++af) {
        for(int hm = 0; hm < ss.hc1_disease_stages; ++hm) {
          if((adj * cpars.hc_art_init_dist(af, time_step)) > 1){
            intermediate.hc_art_scalar = 1;
          }else{
            intermediate.hc_art_scalar = adj * cpars.hc_art_init_dist(af, time_step);
          }
          if(intermediate.hc_art_num > 0){
            intermediate.hc_art_scalar = intermediate.hc_art_scalar;
          }else{
            intermediate.hc_art_scalar =  0.0;
          }
          if(af < ss.hc2_agestart){
            state_next.hc1_art_pop(0, hm, af, g) +=  intermediate.hc_art_scalar * intermediate.hc_art_init(hm, cat, af, g) ;
          }else{
            if(hm < (ss.hc2_disease_stages - 1)){
              state_next.hc2_art_pop(0, hm, af - ss.hc2_agestart, g) +=  intermediate.hc_art_scalar * intermediate.hc_art_init(hm, cat, af, g) ;
            }
          }
          if(af < ss.hc2_agestart){
            state_next.hc1_hiv_pop(hm, cat, af, g) -=  intermediate.hc_art_scalar * intermediate.hc_art_init(hm, cat, af, g) ;
          }else{
            if(hm < (ss.hc2_disease_stages - 1)){
              state_next.hc2_hiv_pop(hm, cat, af - ss.hc2_agestart, g) -=  intermediate.hc_art_scalar * intermediate.hc_art_init(hm, cat, af, g) ;
            }
          }

        } //end ss.hc1_disease_stages
      } // end af
    } // end hTM
  } // end ss.num_genders

}



}// namespace internal

template<HivAgeStratification S, typename real_type>
void run_child_model_simulation(int time_step,
                                const Parameters<real_type> &pars,
                                const State<S, real_type> &state_curr,
                                State<S, real_type> &state_next,
                                internal::IntermediateData<S, real_type> &intermediate) {
  //internal::run_child_ageing(time_step, pars, state_curr, state_next, intermediate);
  //internal::run_child_hiv_infections(time_step, pars, state_curr, state_next, intermediate);
  //internal::run_child_natural_history(time_step, pars, state_curr, state_next, intermediate);
  //internal::run_child_hiv_mort(time_step, pars, state_curr, state_next, intermediate);
 // internal::add_child_grad(time_step, pars, state_curr, state_next, intermediate);
  //this function may need to be broken up, its around 350 lines
  //!!!TODO: also need to fix the looping order for some loops
 // internal::run_child_art_initiation(time_step, pars, state_curr, state_next, intermediate);


}



} // namespace leapfrog

