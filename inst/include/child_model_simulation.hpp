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
  for (int s = 0; s < ss.NS; ++s) {
    //less than 5 because there is a cd4 transition between ages 4 and 5
    for (int a = 1; a < ss.hc2_agestart; ++a) {
      for (int hd = 0; hd < ss.hc1DS; ++hd) {
        for (int cat = 0 ; cat < ss.hcTT; ++cat) {
          state_next.hc1_hiv_pop(hd, cat, a, s) += state_curr.hc1_hiv_pop(hd, cat, a-1, s) * demog.survival_probability(a, s, time_step);

        }
        for (int dur = 0; dur < ss.hTS; ++dur) {
          state_next.hc1_art_pop(dur, hd, a, s) += state_curr.hc1_art_pop(dur, hd, a-1, s) * demog.survival_probability(a, s, time_step);
        }

      }
    }
  }

  for (int s = 0; s < ss.NS; ++s) {
    for (int hd = 0; hd < ss.hc1DS; ++hd) {
      for (int hd_alt = 0; hd_alt < ss.hc2DS; ++hd_alt) {
        for (int cat = 0 ; cat < ss.hcTT; ++cat) {
          state_next.hc2_hiv_pop(hd_alt, cat, 0, s) +=  state_curr.hc1_hiv_pop(hd, cat, ss.hc1_ageend, s) * demog.survival_probability(ss.hc2_agestart, s, time_step) * cpars.hc_cd4_transition(hd_alt, hd);
        }
        for (int dur = 0; dur < ss.hTS; ++dur) {
          state_next.hc2_art_pop(dur, hd_alt, 0, s) += state_curr.hc1_art_pop(dur, hd, ss.hc1_ageend, s) * demog.survival_probability(ss.hc2_agestart, s, time_step) * cpars.hc_cd4_transition(hd_alt, hd);
        }
      }
    }
  }

  for (int s = 0; s < ss.NS; ++s) {
    for (int a = (ss.hc2_agestart + 1); a < pars.options.p_idx_fertility_first; ++a) {
      for (int hd = 0; hd < ss.hc2DS; ++hd) {
        for (int cat = 0 ; cat < ss.hcTT; ++cat) {
          state_next.hc2_hiv_pop(hd, cat, a - ss.hc2_agestart, s) += state_curr.hc2_hiv_pop(hd, cat, a- ss.hc2_agestart-1, s) * demog.survival_probability(a, s, time_step);

        }
        for (int dur = 0; dur < ss.hTS; ++dur) {
          state_next.hc2_art_pop(dur, hd, a - ss.hc2_agestart, s) += state_curr.hc2_art_pop(dur, hd, a- ss.hc2_agestart-1, s) * demog.survival_probability(a, s, time_step);
        }

      }
    }
  }



}

template<HivAgeStratification S, typename real_type>
void calc_hiv_negative_pop(int time_step,
                      const Parameters<real_type> &pars,
                      const State<S, real_type> &state_curr,
                      State<S, real_type> &state_next,
                      IntermediateData<S, real_type> &intermediate) {
  const auto demog = pars.demography;
  const auto cpars = pars.children;
  constexpr auto ss = StateSpace<S>();

  for(int s = 0; s < ss.NS; ++s){
    for(int a = 0; a < ss.hAG; ++a){
      intermediate.p_hiv_neg_pop(a, s) = demog.base_pop(a, s) - state_next.p_hiv_pop(a, s);
    }// end a
  }//end s

}

template<HivAgeStratification S, typename real_type>
void run_perinatal_transmission(int time_step,
                              const Parameters<real_type> &pars,
                              const State<S, real_type> &state_curr,
                              State<S, real_type> &state_next,
                              IntermediateData<S, real_type> &intermediate) {
  constexpr auto ss = StateSpace<S>();
  const auto demog = pars.demography;
  const auto cpars = pars.children;

  intermediate.asfr_sum = 0.0;
  for (int a = 0; a < pars.options.p_fertility_age_groups; ++a) {
    intermediate.asfr_sum += demog.age_specific_fertility_rate(a, time_step);
  } // end a
  //TODO: add in patients reallocated
  //TODO: pull incidence infection mtct into the input object

  for (int hp = 0; hp < ss.hPS; ++hp) {
    intermediate.sumARV += cpars.pmtct(hp,time_step);
  }
  intermediate.number_on_PMTCT = intermediate.sumARV ;// + PATIENTS REALLOCATED

  //Maggie: TODO: make pmtct coverages variables that are automatically converted to % covered, I think that should work?
  //Cap so that it snaps the distribution to one if intermediate.need_pmtcted, with the distribution matching the input numbers
  for (int hp = 0; hp < ss.hPS; ++hp) {
    if (cpars.pmtct_input_is_percent(time_step)) {
      intermediate.pmtct_coverage(hp) = cpars.pmtct(hp, time_step);
    }else{
      //TODO: should this be divided by all pregnant women?
      intermediate.pmtct_coverage(hp) = cpars.pmtct(hp, time_step) / intermediate.sumARV;
    }
  }

  if (cpars.pmtct_input_is_percent(time_step)) { //input is percent
    if ((intermediate.need_pmtctPMTCT(time_step) * intermediate.sumARV) > 0) {
      intermediate.on_pmtct = intermediate.need_pmtctPMTCT(time_step) * intermediate.sumARV ;
    }else{
      intermediate.on_pmtct = 0.0;
    }
  }else{//input is number
    intermediate.on_pmtct = intermediate.sumARV;
  }

  intermediate.need_pmtct = intermediate.need_pmtctPMTCT(time_step) > intermediate.number_on_PMTCT ? intermediate.need_pmtctPMTCT(time_step) : intermediate.number_on_PMTCT;

  intermediate.no_pmtct = (1 - intermediate.sumARV) > 0 ? 1 - intermediate.sumARV : 0;

  //Proportion of pregnant women by CD4 count
  for (int a = 0; a < 35; ++a) {
    intermediate.prop_wlhiv_lt350 += state_curr.h_hiv_adult(4,a,1) + state_curr.h_hiv_adult(5,a,1) + state_curr.h_hiv_adult(6,a,1) ;
    intermediate.num_wlhiv_200to350 += state_curr.h_hiv_adult(3,a,1) + state_curr.h_hiv_adult(2,a,1) ;
    intermediate.num_wlhiv_gte350 += state_curr.h_hiv_adult(0,a,1) + state_curr.h_hiv_adult(1,a,1) ;
  }

  intermediate.num_wlhiv = intermediate.num_wlhiv_200to350 + intermediate.num_wlhiv_gte350 + intermediate.prop_wlhiv_lt350;


  if (intermediate.num_wlhiv >0) {
    intermediate.prop_wlhiv_lt200 = intermediate.prop_wlhiv_lt350/ intermediate.num_wlhiv;
    intermediate.prop_wlhiv_200to350 = intermediate.num_wlhiv_200to350 / intermediate.num_wlhiv;
    intermediate.prop_wlhiv_gte350 = intermediate.num_wlhiv_gte350 / intermediate.num_wlhiv;
  }else{
    intermediate.prop_wlhiv_lt200 = 0;
    intermediate.prop_wlhiv_200to350 = 1;
    intermediate.prop_wlhiv_gte350 = 0;
  }

  intermediate.prop_wlhiv_lt350 = intermediate.prop_wlhiv_lt200 + intermediate.prop_wlhiv_200to350;
  //adjust option A and option B
  //Option A and B were only authorized for women with greater than 350 CD4, so if the percentage of women
  //on option A/B > the proportion of women in this cd4 category, we assume that some must have a cd4 less than 350
  //option AB will be less effective for these women so we adjust for that
  if ((cpars.pmtct(0,time_step) + cpars.pmtct(1,time_step)) > intermediate.prop_wlhiv_gte350) {
    if (intermediate.prop_wlhiv_gte350 > 0) {
      intermediate.excessratio = ((cpars.pmtct(0,time_step) + cpars.pmtct(1,time_step)) / intermediate.prop_wlhiv_gte350) - 1;
    }else{
      intermediate.excessratio = 0;
    }
    intermediate.optA_transmission_rate = cpars.pmtct_transmission_rate(0,1,0) * (1 + intermediate.excessratio);
    intermediate.optB_transmission_rate = cpars.pmtct_transmission_rate(0,2,0) * (1 + intermediate.excessratio);
  }

  else{
    intermediate.excessratio = 0.0;
    intermediate.optA_transmission_rate = cpars.pmtct_transmission_rate(0,1,0) * (1 + intermediate.excessratio);
    intermediate.optB_transmission_rate = cpars.pmtct_transmission_rate(0,2,0) * (1 + intermediate.excessratio);
  }


  ///////////////////////////////////
  //Calculate transmission rate
  ///////////////////////////////////
  intermediate.retained_on_ART = cpars.pmtct(4,time_step) * cpars.pmtct_dropout(4,time_step);
  intermediate.retained_started_ART = cpars.pmtct(5,time_step) * cpars.pmtct_dropout(5,time_step);
  //Transmission among women on treatment
  intermediate.perinatal_transmission_rate = cpars.pmtct(2,time_step) * cpars.pmtct_transmission_rate(0,2,0) + cpars.pmtct(3,time_step) * cpars.pmtct_transmission_rate(0,3,0) + cpars.pmtct(0,time_step) * intermediate.optA_transmission_rate + cpars.pmtct(1,time_step) * intermediate.optB_transmission_rate + intermediate.retained_on_ART * cpars.pmtct_transmission_rate(0,4,0) + intermediate.retained_started_ART * cpars.pmtct_transmission_rate(0,5,0)+ cpars.pmtct(6,time_step) * cpars.pmtct_transmission_rate(0,6,0);

  intermediate.receiving_pmtct = cpars.pmtct(0,time_step) + cpars.pmtct(1,time_step) + cpars.pmtct(2,time_step) + cpars.pmtct(3,time_step) + intermediate.retained_on_ART + intermediate.retained_started_ART + cpars.pmtct(6,time_step);
  intermediate.no_pmtct = 1 - intermediate.receiving_pmtct;

  //Transmission among women not on treatment
  if (intermediate.num_wlhiv > 0) {
    intermediate.perinatal_transmission_rate = intermediate.perinatal_transmission_rate + intermediate.no_pmtct * (intermediate.prop_wlhiv_lt200 * cpars.vertical_transmission_rate(4) + intermediate.prop_wlhiv_200to350 * cpars.vertical_transmission_rate(2) + intermediate.prop_wlhiv_gte350 * cpars.vertical_transmission_rate(0));
  }else{
    intermediate.perinatal_transmission_rate = intermediate.perinatal_transmission_rate;
  }
  intermediate.perinatal_transmission_rate_bf_calc = intermediate.perinatal_transmission_rate;

  //Transmission due to incident infections
  for (int a = 0; a < 35; ++a) {
    intermediate.age_weighted_hivneg += demog.age_specific_fertility_rate(a, time_step) / intermediate.asfr_sum  * intermediate.p_hiv_neg_pop(a + 15,1) ; //HIV negative 15-49 women weighted for ASFR
    intermediate.age_weighted_infections +=  demog.age_specific_fertility_rate(a, time_step) / intermediate.asfr_sum  * state_curr.infections(a + 15,1) ; //newly infected 15-49 women, weighted for ASFR
  }

  intermediate.incidence_rate_wlhiv = intermediate.age_weighted_infections / intermediate.age_weighted_hivneg;
  intermediate.perinatal_transmission_from_incidence = intermediate.incidence_rate_wlhiv * (9/12) * (intermediate.births_sum - intermediate.need_pmtct) *0;  //.181;

  if (intermediate.need_pmtct > 0) {
    intermediate.perinatal_transmission_rate = intermediate.perinatal_transmission_rate + intermediate.perinatal_transmission_from_incidence / intermediate.need_pmtct;
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

  for (int s = 0; s < ss.NS; ++s) {
    // Run only first 5 age groups in total population 0, 1, 2, 3, 4
    for (int a = 0; a < ss.hc2_agestart; ++a) {
      if (cpars.hc_nosocomial(time_step) > 0) {
        state_next.p_infections(a, s) = cpars.hc_nosocomial(time_step) / (5.0 * ss.NS); //5.0 is used because we want to evenly distribute across the 5 age groups in 0-4
        state_next.p_hiv_pop(a, s) += state_next.p_infections(a, s);

        for (int hd = 0; hd < ss.hc1DS; ++hd) {
          // putting them all in perinatal hcTT to match spec nosocomial
          if (cpars.hc1_cd4_dist(hd) > 0) {
            state_next.hc1_hiv_pop(hd, 0, a, s) += state_next.p_infections(a, s) * cpars.hc1_cd4_dist(hd);
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


 for (int s = 0; s < ss.NS; ++s) {
   for (int hd = 0; hd < ss.hc1DS; ++hd) {
     for (int a = 0; a < ss.hc2_agestart; ++a) {
       for (int cat = 0; cat < ss.hcTT; ++cat) {
         intermediate.hc_posthivmort(hd, cat, a, s) += state_next.hc1_hiv_pop(hd, cat, a, s) - (1.0 - cpars.ctx_effect * cpars.ctx_val(time_step)) * state_curr.hc1_hiv_pop(hd, cat, a, s) * cpars.hc1_cd4_mort(hd, cat, a);
          state_next.hc1_noart_aids_deaths(hd, cat, a, s) += (1 - cpars.ctx_effect * cpars.ctx_val(time_step)) * state_curr.hc1_hiv_pop(hd, cat, a, s) * cpars.hc1_cd4_mort(hd, cat, a)  ;
       }
     }
   }
 }


 for (int s = 0; s < ss.NS; ++s) {
   for (int hd = 0; hd < ss.hc2DS; ++hd) {
     for (int a = ss.hc2_agestart; a < pars.options.p_idx_fertility_first; ++a) {
       for (int cat = 0; cat < ss.hcTT; ++cat) {
         intermediate.hc_posthivmort(hd, cat, a, s) += state_next.hc2_hiv_pop(hd, cat, a - ss.hc2_agestart, s) - (1 - cpars.ctx_effect * cpars.ctx_val(time_step)) * state_curr.hc2_hiv_pop(hd, cat, a - ss.hc2_agestart, s) * cpars.hc2_cd4_mort(hd, cat, a - ss.hc2_agestart);
         state_next.hc2_noart_aids_deaths(hd, cat, a, s) +=  (1 - cpars.ctx_effect * cpars.ctx_val(time_step)) * state_curr.hc2_hiv_pop(hd, cat, a, s) * cpars.hc2_cd4_mort(hd, cat, a - 5); // output hiv deaths, aggregated across transmission category
       }
     }
   }
 }

 //progress through CD4 categories
 for (int s = 0; s < ss.NS; ++s) {
   for (int hd = 1; hd < ss.hc1DS; ++hd) {
     for (int a = 0; a < ss.hc2_agestart; ++a) {
       for (int cat = 0; cat < ss.hcTT; ++cat) {
         intermediate.hc_grad(hd - 1, cat, a, s) -=  (intermediate.hc_posthivmort(hd - 1, cat, a, s) * cpars.hc1_cd4_prog(hd - 1) + state_next.hc1_hiv_pop(hd - 1, cat, a, s) * cpars.hc1_cd4_prog(hd - 1)) / 2; //moving to next cd4 category
         intermediate.hc_grad(hd, cat, a, s) += (intermediate.hc_posthivmort(hd - 1, cat, a, s) * cpars.hc1_cd4_prog(hd - 1) + state_next.hc1_hiv_pop(hd - 1, cat, a, s) * cpars.hc1_cd4_prog(hd - 1)) / 2; //moving into this cd4 category
       }
     }
   }
 }

 //progress through CD4 categories
 for (int s = 0; s < ss.NS; ++s) {
   for (int hd = 1; hd < ss.hc2DS; ++hd) {
     for (int a = ss.hc2_agestart; a < pars.options.p_idx_fertility_first; ++a) {
       for (int cat = 0; cat < ss.hcTT; ++cat) {
         intermediate.hc_grad(hd - 1, cat, a, s) -= (intermediate.hc_posthivmort(hd - 1, cat, a, s) * cpars.hc2_cd4_prog(hd - 1) + state_next.hc2_hiv_pop(hd - 1, cat, a - ss.hc2_agestart, s) * cpars.hc2_cd4_prog(hd - 1)) / 2; //moving to next cd4 category
         intermediate.hc_grad(hd, cat, a, s) += (intermediate.hc_posthivmort(hd - 1, cat, a, s) * cpars.hc2_cd4_prog(hd - 1) + state_next.hc2_hiv_pop(hd - 1, cat, a - ss.hc2_agestart, s) * cpars.hc2_cd4_prog(hd - 1)) / 2; //moving into this cd4 category
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


  for (int s = 0; s < ss.NS; ++s) {
    for (int hd = 0; hd < ss.hc2DS; ++hd) {
      for (int a = 0; a < pars.options.p_idx_fertility_first; ++a) {
        for (int cat = 0; cat < ss.hcTT; ++cat) {
          if(a < ss.hc2_agestart){
            intermediate.hc_grad(hd, cat, a, s) -= (1 - cpars.ctx_effect * cpars.ctx_val(time_step)) * state_next.hc1_hiv_pop(hd, cat, a, s) * cpars.hc1_cd4_mort(hd, cat, a)  ;
          }else{
            intermediate.hc_grad(hd, cat, a, s) -= (1 - cpars.ctx_effect * cpars.ctx_val(time_step)) * state_next.hc2_hiv_pop(hd, cat, a - ss.hc2_agestart, s) * cpars.hc2_cd4_mort(hd, cat, a - ss.hc2_agestart)  ;

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
  for (int s = 0; s < ss.NS; ++s) {
    for (int hd = 0; hd < ss.hc1DS; ++hd) {
      for (int a = 0; a < pars.options.p_idx_fertility_first; ++a) {
        for (int cat = 0; cat < ss.hcTT; ++cat) {
          if(a < ss.hc2_agestart){
            state_next.hc1_hiv_pop(hd, cat, a, s) +=intermediate.hc_grad(hd, cat, a, s);
          }else{
          state_next.hc2_hiv_pop(hd, cat, a- ss.hc2_agestart, s) +=intermediate.hc_grad(hd, cat, a, s);
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
  for(int s = 0; s <ss.NS; ++s)  {
    for(int cat = 0; cat < ss.hcTT; ++cat) {
      for(int a = 0; a < pars.options.p_idx_fertility_first; ++a) {
        for(int hd = 0; hd < ss.hc1DS; ++hd) {
          if(a < cpars.hc_art_elig_age(time_step)){
            if(a < ss.hc2_agestart){
              intermediate.hc_art_need(hd, cat, a, s) += state_next.hc1_hiv_pop(hd, cat, a, s);
            }else{
            if(hd < ss.hc2DS){
              intermediate.hc_art_need(hd, cat, a, s) += state_next.hc2_hiv_pop(hd, cat, a - ss.hc2_agestart, s);
            }
            }
          }
        } // end ss.hc1DS
      } // end a
    } // end hcTT
  } // end ss.NS

  //all children under a certain CD4 eligible for ART
  for(int s = 0; s <ss.NS; ++s)  {
    for(int cat = 0; cat < ss.hcTT; ++cat) {
      for(int a = cpars.hc_art_elig_age(time_step); a < pars.options.p_idx_fertility_first; ++a) {
        for(int hd = 0; hd < ss.hc1DS; ++hd) {
          if(hd > (cpars.hc_art_elig_cd4(a, time_step) - 2)){ //-2 to account for the zero-based indexing and basically as >=
            if(a < ss.hc2_agestart){
              intermediate.hc_art_need(hd, cat, a, s) += state_next.hc1_hiv_pop(hd, cat, a, s);
            }else{
              if(hd < ss.hc2DS){
                intermediate.hc_art_need(hd, cat, a, s) += state_next.hc2_hiv_pop(hd, cat, a - ss.hc2_agestart, s);
              }
            }
          }
        } // end ss.hc1DS
      } // end a
    } // end hcTT
  } // end ss.NS

  for(int s = 0; s <ss.NS; ++s) {
    for(int a = 0; a < pars.options.p_idx_fertility_first; ++a) {
      for(int hd = 0; hd < ss.hc1DS; ++hd) {
        for(int cat = 0; cat < ss.hcTT; ++cat) {
          state_next.hc_art_num += intermediate.hc_art_need(hd, cat, a, s);
        } // end hcTT
      } // end ss.hc1DS
    } // end a
  } // end ss.NS

  //how many should initialize ART
  for(int s = 0; s <ss.NS; ++s) {
    for(int a = 0; a < pars.options.p_idx_fertility_first; ++a) {
      for(int hd = 0; hd < ss.hc1DS; ++hd) {
        for(int cat = 0; cat < ss.hcTT; ++cat) {
          if(cpars.hc_art_val(time_step) > 0){
            intermediate.hc_art_init(hd, cat, a, s) += intermediate.hc_art_need(hd, cat, a, s);
            for(int dur = 0; dur < ss.hTS; ++dur) {
              if(intermediate.hc_art_init(hd, cat, a, s) < 0){
                intermediate.hc_art_init(hd, cat, a, s) = 0.0;
              }else{
                intermediate.hc_art_init(hd, cat, a, s) = intermediate.hc_art_init(hd, cat, a, s);
              }
            }// end ss.hTS
          }
        }// end hcTT
      }// end ss.hc1DS
    }// end a
  }// end ss.NS


  for(int s = 0; s <ss.NS; ++s) {
    for(int a = 0; a < pars.options.p_idx_fertility_first; ++a) {
      for(int hd = 0; hd < ss.hc1DS; ++hd) {
        for(int cat = 0; cat < ss.hcTT; ++cat) {
          intermediate.hc_art_init_total += intermediate.hc_art_init(hd, cat, a, s);
        }// end hcTT
      }// end ss.hc1DS
    }// end a
  }// end ss.NS

  //!!! TODO: fix order of for loop
  for(int s = 0; s <ss.NS; ++s) {
    for(int hd = 0; hd < ss.hc1DS; ++hd) {
      for(int a = 0; a < pars.options.p_idx_fertility_first; ++a) {
        intermediate.hc_death_rate = 0.0;
        intermediate.hc_art_grad(0, hd, a, s) = 0.0;

        if(a < ss.hc2_agestart){
          if(state_next.hc1_art_pop(0, hd, a, s) >0){
            intermediate.hc_death_rate =  cpars.hc_art_mort_rr(0, a, time_step) * 0.5 * (cpars.hc1_art_mort(hd, 0, a) + cpars.hc1_art_mort(hd, 1, a)) ;
            state_next.hc1_art_aids_deaths(0,hd, a, s) =  intermediate.hc_death_rate * state_next.hc1_art_pop(0, hd, a, s)  ;
            intermediate.hc_art_grad(0,hd, a, s) -= state_next.hc1_art_aids_deaths(0,hd, a, s) ;
            state_next.hc1_art_pop(0, hd,  a, s) += intermediate.hc_art_grad(0, hd, a, s) ;
          }
        }else{
          if(hd < ss.hc2DS){
            if(state_next.hc2_art_pop(0, hd, a-ss.hc2_agestart, s) >0){
            intermediate.hc_death_rate =  cpars.hc_art_mort_rr(0, a, time_step) * 0.5 * (cpars.hc2_art_mort(hd, 0, a-ss.hc2_agestart) + cpars.hc2_art_mort(hd, 1, a-ss.hc2_agestart));
            state_next.hc2_art_aids_deaths(0,hd, a-ss.hc2_agestart, s) =  intermediate.hc_death_rate * state_next.hc2_art_pop(0, hd, a-ss.hc2_agestart, s)  ;
            intermediate.hc_art_grad(0, hd, a, s) -= state_next.hc2_art_aids_deaths(0,hd, a-ss.hc2_agestart, s) ;
            state_next.hc2_art_pop(0, hd,  a-ss.hc2_agestart, s) += intermediate.hc_art_grad(0, hd, a, s) ;
          }
         }
        }
      }// end a
    }// end ss.hc1DS
  }// end ss.NS


  for(int s = 0; s <ss.NS; ++s) {
    for(int hd = 0; hd < ss.hc1DS; ++hd) {
      for(int a = 0; a < pars.options.p_idx_fertility_first; ++a) {
        intermediate.hc_death_rate = 0;
        intermediate.hc_art_grad(2, hd, a, s) = 0.0;

        if(a < ss.hc2_agestart){
          intermediate.hc_death_rate = cpars.hc_art_mort_rr(2, a, time_step) * cpars.hc1_art_mort(hd, 2, a);
          state_next.hc1_art_aids_deaths(2,hd, a, s) =  state_next.hc1_art_pop(2, hd, a, s) * intermediate.hc_death_rate;
          intermediate.hc_art_grad(2,hd, a, s) -= state_next.hc1_art_aids_deaths(2,hd, a, s) ;
          state_next.hc1_art_pop(2, hd,  a, s) += intermediate.hc_art_grad(2, hd, a, s) ;
        }else{
          if(hd < (ss.hc2DS)){
            intermediate.hc_death_rate = cpars.hc_art_mort_rr(2, a, time_step) * cpars.hc2_art_mort(hd, 2, a-ss.hc2_agestart);
            state_next.hc2_art_aids_deaths(2,hd, a-ss.hc2_agestart, s) =  state_next.hc2_art_pop(2, hd, a-ss.hc2_agestart, s) * intermediate.hc_death_rate;
            intermediate.hc_art_grad(2,hd, a, s) -= state_next.hc2_art_aids_deaths(2,hd, a-ss.hc2_agestart, s) ;
            state_next.hc2_art_pop(2, hd,  a-ss.hc2_agestart, s) += intermediate.hc_art_grad(2, hd, a, s) ;
          }
        }
      }// end a
    }// end ss.hc1DS
  }// end ss.NS

  //Progress ART to the correct time on ART
  //!!! TODO: fix order of for loop
  for(int hd = 0; hd < ss.hc1DS; ++hd) {
    for(int a = 0; a < pars.options.p_idx_fertility_first; ++a) {
      for(int s = 0; s <ss.NS; ++s) {
        if(a < ss.hc2_agestart){
          if(state_next.hc1_art_pop(0, hd, a, s) > 0){
            state_next.hc1_art_pop(1, hd, a, s) += state_next.hc1_art_pop(0, hd, a, s);
            state_next.hc1_art_pop(0, hd, a, s) -= state_next.hc1_art_pop(0, hd, a, s);
          }
        }else{
        //  if(state_next.hc2_art_pop(0, hd, a-ss.hc2_agestart, s) > 0){
            if(hd < (ss.hc2DS)){
              state_next.hc2_art_pop(1, hd, a-ss.hc2_agestart, s) += state_next.hc2_art_pop(0, hd, a-ss.hc2_agestart, s);
              state_next.hc2_art_pop(0, hd, a-ss.hc2_agestart, s) -= state_next.hc2_art_pop(0, hd, a-ss.hc2_agestart, s);
            }
        //  }
        }
      }//end ss.NS
    }// end a
  }// end ss.hc1DS


  if ( !cpars.hc_art_isperc(time_step) && !cpars.hc_art_isperc(time_step-1) ){ // both numbers
    //Remove how many that are already on ART
    state_next.hc_art_num =  (cpars.hc_art_val(time_step) + cpars.hc_art_val(time_step-1)) / 2 ;
    for(int s = 0; s <ss.NS; ++s) {
      for(int a = 0; a < pars.options.p_idx_fertility_first; ++a) {
        for(int hd = 0; hd < ss.hc1DS; ++hd) {
          for(int dur = 0; dur < ss.hTS; ++dur) {
            if(a < ss.hc2_agestart){
              state_next.hc_art_num -= state_next.hc1_art_pop(dur, hd, a, s)  ;
            }else{
              if(hd < (ss.hc2DS)){
                state_next.hc_art_num -= state_next.hc2_art_pop(dur, hd, a-ss.hc2_agestart, s)   ;
              }
            }
          }// end ss.hTS
        }// end ss.hc1DS
      }// end a
    }// end ss.NS

    if(intermediate.hc_art_init_total < state_next.hc_art_num){
      state_next.hc_art_num = intermediate.hc_art_init_total;
    }
    if(state_next.hc_art_num < 0){
      state_next.hc_art_num =  0;
    }

  } else if (cpars.hc_art_isperc(time_step) && cpars.hc_art_isperc(time_step-1)){ // both percentages
    for(int s = 0; s <ss.NS; ++s) {
      for(int a = 0; a < pars.options.p_idx_fertility_first; ++a) {
        for(int hd = 0; hd < ss.hc1DS; ++hd) {
          for(int dur = 0; dur < ss.hTS; ++dur) {
            if(a < ss.hc2_agestart){
              state_next.hc_art_num += state_next.hc1_art_pop(dur, hd, a, s)  ;
              state_next.hc_art_num += state_next.hc1_art_aids_deaths(dur,hd, a, s) ;
            }else{
              if(hd < (ss.hc2DS )){
                state_next.hc_art_num += state_next.hc2_art_pop(dur, hd, a-ss.hc2_agestart, s)  ;
                state_next.hc_art_num += state_next.hc2_art_aids_deaths(dur,hd, a-ss.hc2_agestart, s) ;
              }
            }
          } //end ss.hTS
        } // end ss.hC1_disease_stages
      } // end a
    } // end ss.NS

    state_next.hc_art_num =  state_next.hc_art_num * (cpars.hc_art_val(time_step) + cpars.hc_art_val(time_step-1)) / 2 ;

    //Remove how many that are already on ART
    for(int s = 0; s <ss.NS; ++s) {
      for(int a = 0; a < pars.options.p_idx_fertility_first; ++a) {
        for(int hd = 0; hd < ss.hc1DS; ++hd) {
          for(int dur = 0; dur < ss.hTS; ++dur) {
            if(a < ss.hc2_agestart){
              state_next.hc_art_num -= state_next.hc1_art_pop(dur, hd, a, s)  ;
            }else{
              if(hd < (ss.hc2DS)){
                state_next.hc_art_num -= state_next.hc2_art_pop(dur, hd, a-ss.hc2_agestart, s)  ;
              }
            }
          } // end ss.hTS
        } // end ss.hC1_disease_stages
      } // end a
    } // end ss.NS
    if(intermediate.hc_art_init_total < state_next.hc_art_num){
      state_next.hc_art_num = intermediate.hc_art_init_total;
    }



  } else if (cpars.hc_art_isperc(time_step) && !cpars.hc_art_isperc(time_step-1)){ // num to percentage
    //Remove how many that are already on ART
    for(int s = 0; s <ss.NS; ++s) {
      for(int a = 0; a < pars.options.p_idx_fertility_first; ++a) {
        for(int hd = 0; hd < ss.hc1DS; ++hd) {
          for(int dur = 0; dur < ss.hTS; ++dur) {
            if(a < ss.hc2_agestart){
              state_next.hc_art_num += state_next.hc1_art_pop(dur, hd, a, s) +  state_next.hc1_art_aids_deaths(dur,hd, a, s)  ;
            }else{
              if(hd < (ss.hc2DS)){
                state_next.hc_art_num += state_next.hc2_art_pop(dur, hd, a-ss.hc2_agestart, s) +  state_next.hc2_art_aids_deaths(dur,hd, a-ss.hc2_agestart, s)  ;
              }
            }
          } // end ss.hTS
        } // end ss.hc1DS
      } //end a
    } //end ss.NS
    state_next.hc_art_num = (cpars.hc_art_val(time_step-1) + (state_next.hc_art_num * cpars.hc_art_val(time_step))) / 2 ;

    for(int s = 0; s <ss.NS; ++s) {
      for(int a = 0; a < pars.options.p_idx_fertility_first; ++a) {
        for(int hd = 0; hd < ss.hc1DS; ++hd) {
          for(int dur = 0; dur < ss.hTS; ++dur) {
            if(a < ss.hc2_agestart){
              state_next.hc_art_num -= state_next.hc1_art_pop(dur, hd, a, s);
            }else{
              if(hd < (ss.hc2DS )){
                state_next.hc_art_num -= state_next.hc2_art_pop(dur, hd, a-ss.hc2_agestart, s);
              }
            }
          } // end ss.hTS
        } // end ss.hc1DS
      } // end a
    } // end ss.NS

    if(state_next.hc_art_num < 0){
      state_next.hc_art_num = 0;
    }
    if(intermediate.hc_art_init_total < state_next.hc_art_num){
      state_next.hc_art_num = intermediate.hc_art_init_total;
    }

  } else if (cpars.hc_art_isperc(time_step-1) && !cpars.hc_art_isperc(time_step)){ //percentage to num


    for(int s = 0; s <ss.NS; ++s) {
      for(int a = 0; a < pars.options.p_idx_fertility_first; ++a) {
        for(int hd = 0; hd < ss.hc1DS; ++hd) {
          for(int dur = 0; dur < ss.hTS; ++dur) {
            if(a < ss.hc2_agestart){
              state_next.hc_art_num -= state_next.hc1_art_pop(dur, hd, a, s)  ;
            }else{
              if(hd < (ss.hc2DS)){
                state_next.hc_art_num -= state_next.hc2_art_pop(dur, hd, a-ss.hc2_agestart, s)   ;
              }
            }
          } // end ss.hTS
        } //end ss.hc1DS
      } // end a
    } //end ss.NS


   state_next.hc_art_num = (state_curr.hc_art_num + cpars.hc_art_val(time_step)) / 2 ;

    //Remove how many that are already on ART
    for(int s = 0; s <ss.NS; ++s) {
      for(int a = 0; a < pars.options.p_idx_fertility_first; ++a) {
        for(int hd = 0; hd < ss.hc1DS; ++hd) {
          for(int dur = 0; dur < ss.hTS; ++dur) {
            if(a < ss.hc2_agestart){
              state_next.hc_art_num -= state_next.hc1_art_pop(dur, hd, a, s);
            }else{
              if(hd < (ss.hc2DS)){
                state_next.hc_art_num -= state_next.hc2_art_pop(dur, hd, a-ss.hc2_agestart, s);
              }
            }
          } //end ss.hTS
        } //end ss.hc1DS
      } //end a
    } //end ss.NS

    if(state_next.hc_art_num < 0){
      state_next.hc_art_num = 0;
    }
    if(intermediate.hc_art_init_total < state_next.hc_art_num){
      state_next.hc_art_num = intermediate.hc_art_init_total;
    }
  }



  for(int s = 0; s <ss.NS; ++s) {
    for(int a = 0; a < pars.options.p_idx_fertility_first; ++a) {
      for(int hd = 0; hd < ss.hc1DS; ++hd) {
        for(int cat = 0; cat < ss.hcTT; ++cat) {
          intermediate.hc_initByAge +=  intermediate.hc_art_init(hd, cat, a, s) * cpars.hc_art_init_dist(a, time_step);
        } //end hcTT
      } // end ss.hc1DS
    } //end a
  } // end ss.NS

  if(intermediate.hc_initByAge == 0.0){
    intermediate.hc_adj = 1.0 ;
  }else{
    intermediate.hc_adj = state_next.hc_art_num / intermediate.hc_initByAge ;
  }
  for(int s = 0; s <ss.NS; ++s) {
    for(int cat = 0; cat < ss.hcTT; ++cat) {
      for(int a = 0; a < pars.options.p_idx_fertility_first; ++a) {
        for(int hd = 0; hd < ss.hc1DS; ++hd) {
          if((intermediate.hc_adj * cpars.hc_art_init_dist(a, time_step)) > 1.0){
            intermediate.hc_art_scalar = 1.0;
          }else{
            intermediate.hc_art_scalar = intermediate.hc_adj * cpars.hc_art_init_dist(a, time_step);
          }
          if(state_next.hc_art_num > 0.0){
            intermediate.hc_art_scalar = intermediate.hc_art_scalar;
          }else{
            intermediate.hc_art_scalar =  0.0;
          }
          if(a < ss.hc2_agestart){
            state_next.hc1_art_pop(0, hd, a, s) +=  intermediate.hc_art_scalar * intermediate.hc_art_init(hd, cat, a, s) ;
          }else{
            if(hd < (ss.hc2DS)){
              state_next.hc2_art_pop(0, hd, a - ss.hc2_agestart, s) +=  intermediate.hc_art_scalar * intermediate.hc_art_init(hd, cat, a, s) ;
            }
          }
          if(a < ss.hc2_agestart){
            state_next.hc1_hiv_pop(hd, cat, a, s) -=  intermediate.hc_art_scalar * intermediate.hc_art_init(hd, cat, a, s) ;
          }else{
            if(hd < (ss.hc2DS )){
              state_next.hc2_hiv_pop(hd, cat, a - ss.hc2_agestart, s) -=  intermediate.hc_art_scalar * intermediate.hc_art_init(hd, cat, a, s) ;
            }
          }

        } //end ss.hc1DS
      } // end a
    } // end hcTT
  } // end ss.NS
}


template<HivAgeStratification S, typename real_type>
void run_child_art_mortality(int time_step,
                           const Parameters<real_type> &pars,
                           const State<S, real_type> &state_curr,
                           State<S, real_type> &state_next,
                           IntermediateData<S, real_type> &intermediate) {
  constexpr auto ss = StateSpace<S>();
  const auto cpars = pars.children;
  intermediate.hc_art_grad.setZero();

  for(int hd = 0; hd < ss.hc1DS; ++hd) {
    for(int a = 0; a < pars.options.p_idx_fertility_first; ++a) {
      for(int s = 0; s <ss.NS; ++s) {
        intermediate.hc_death_rate = 0.0;
        intermediate.hc_art_grad(0, hd, a, s) = 0.0;
        if(a < ss.hc2_agestart){
          if(state_next.hc1_art_pop(0, hd, a, s) > 0){
            intermediate.hc_death_rate = cpars.hc_art_mort_rr(0, a, time_step) * 0.5 * (cpars.hc1_art_mort(hd, 0, a) + cpars.hc1_art_mort(hd, 1, a));
            state_next.hc1_art_aids_deaths(0, hd, a, s) = intermediate.hc_death_rate * state_next.hc1_art_pop(0, hd, a, s);
            intermediate.hc_art_grad(0, hd, a, s) -= state_next.hc1_art_aids_deaths(0, hd, a, s);
            state_next.hc1_art_pop(0, hd, a, s) += intermediate.hc_art_grad(0, hd, a, s);
          }
        }else{
          if(hd < ss.hc2DS){
            intermediate.hc_death_rate = cpars.hc_art_mort_rr(0, a, time_step) * 0.5 * (cpars.hc2_art_mort(hd, 0, a-ss.hc2_agestart) + cpars.hc2_art_mort(hd, 1, a-ss.hc2_agestart));
            state_next.hc2_art_aids_deaths(0, hd, a-ss.hc2_agestart, s)  = intermediate.hc_death_rate * state_next.hc2_art_pop(0, hd, a-ss.hc2_agestart, s);
            intermediate.hc_art_grad(0, hd, a, s) -= state_next.hc2_art_aids_deaths(0, hd, a-ss.hc2_agestart, s);
            state_next.hc2_art_pop(0, hd, a-ss.hc2_agestart, s) += intermediate.hc_art_grad(0, hd, a, s);
          }
        }
      }//end ss.NS
    }// end a
  }// end ss.hc1DS

  for(int hd = 0; hd < ss.hc1DS; ++hd) {
    for(int a = 0; a < pars.options.p_idx_fertility_first; ++a) {
      for(int s = 0; s <ss.NS; ++s) {
        if(a < ss.hc2_agestart){
          if(state_next.hc1_art_pop(1, hd, a, s) > 0){
            state_next.hc1_art_pop(2, hd, a, s) += state_next.hc1_art_pop(1, hd, a, s);
            state_next.hc1_art_pop(1, hd, a, s) -= state_next.hc1_art_pop(1, hd, a, s);
          }
        }else{
          if(state_next.hc2_art_pop(1, hd, a-ss.hc2_agestart, s) > 0){
            if(hd < (ss.hc2DS)){
              state_next.hc2_art_pop(2, hd, a-ss.hc2_agestart, s) += state_next.hc2_art_pop(1, hd, a-ss.hc2_agestart, s);
              state_next.hc2_art_pop(1, hd, a-ss.hc2_agestart, s) -= state_next.hc2_art_pop(1, hd, a-ss.hc2_agestart, s);
            }
          }
        }
      }//end ss.NS
    }// end a
  }// end ss.hc1DS
}

template<HivAgeStratification S, typename real_type>
void run_wlhiv_births(int time_step,
                             const Parameters<real_type> &pars,
                             const State<S, real_type> &state_curr,
                             State<S, real_type> &state_next,
                             IntermediateData<S, real_type> &intermediate) {
  constexpr auto ss = StateSpace<S>();
  const auto demog = pars.demography;
  const auto cpars = pars.children;


  for (int a = 0; a < pars.options.p_fertility_age_groups; ++a) {
    intermediate.asfr_sum += demog.age_specific_fertility_rate(a, time_step);
  } // end a

  intermediate.births_sum = state_next.births;


  for (int a = 0; a < pars.options.p_fertility_age_groups; ++a) {
    intermediate.nHIVcurr = 0.0;
    intermediate.nHIVlast = 0.0;
    intermediate.df = 0.0;

    for (int hd = 0; hd < ss.hDS; ++hd) {
      intermediate.nHIVcurr += state_next.h_hiv_adult(hd, a, 1);
      intermediate.nHIVlast += state_curr.h_hiv_adult(hd, a, 1);
      for (int ht = 0; ht < ss.hTS; ++ht) {
        intermediate.nHIVcurr += state_next.h_art_adult(ht, hd, a, 1);
        intermediate.nHIVlast += state_curr.h_art_adult(ht, hd, a, 1);
      } //end hTS
    } //end hDS

   // intermediate.totpop = intermediate.nHIVcurr + state_next.hivnpop1(a + 15,1);
    intermediate.prev = intermediate.nHIVcurr / state_next.p_total_pop(a + 15,1);


    for (int hd = 0; hd < ss.hDS; ++hd) {
      intermediate.df += cpars.local_adj_factor * cpars.fert_mult_by_age(a) * cpars.fert_mult_offart(hd) * ((state_next.h_hiv_adult(hd, a, 1) + state_curr.h_hiv_adult(hd, a, 1)) / 2);
      //women on ART less than 6 months use the off art fertility multiplier
      intermediate.df += cpars.local_adj_factor * cpars.fert_mult_by_age(a) * cpars.fert_mult_offart(hd) * ((state_next.h_art_adult(0, hd, a, 1) + state_curr.h_art_adult(0, hd, a, 1)) / 2);
      for (int ht = 1; ht < ss.hTS; ++ht) {
        intermediate.df += cpars.local_adj_factor * cpars.fert_mult_onart(a) * ((state_next.h_art_adult(ht, hd, a, 1) + state_curr.h_art_adult(ht, hd, a, 1)) / 2);
      } //end hTS
    } // end hDS


    if (intermediate.nHIVcurr > 0) {
      intermediate.df = intermediate.df / ((intermediate.nHIVcurr + intermediate.nHIVlast) / 2);
    }else{
      intermediate.df = 1;
    }


    intermediate.birthsCurrAge = (intermediate.nHIVcurr + intermediate.nHIVlast) / 2 * cpars.total_fertility_rate(time_step) * intermediate.df / (intermediate.df * intermediate.prev + 1 - intermediate.prev) *  demog.age_specific_fertility_rate(a, time_step) / intermediate.asfr_sum ;
    intermediate.birthsHE += intermediate.birthsCurrAge;
    if (a < 9) {
      intermediate.births_HE_15_24 += intermediate.birthsCurrAge;
    }

  } // end a

  state_next.hiv_births = intermediate.birthsHE;

}


}// namespace internal

template<HivAgeStratification S, typename real_type>
void run_child_model_simulation(int time_step,
                                const Parameters<real_type> &pars,
                                const State<S, real_type> &state_curr,
                                State<S, real_type> &state_next,
                                internal::IntermediateData<S, real_type> &intermediate) {
 internal::run_child_ageing(time_step, pars, state_curr, state_next, intermediate);
 internal::run_child_hiv_infections(time_step, pars, state_curr, state_next, intermediate);
 internal::run_child_natural_history(time_step, pars, state_curr, state_next, intermediate);
 internal::run_child_hiv_mort(time_step, pars, state_curr, state_next, intermediate);
 internal::add_child_grad(time_step, pars, state_curr, state_next, intermediate);
 // this function may intermediate.need_pmtct to be broken up, its around 350 lines
 // !!!TODO: also intermediate.need_pmtct to fix the looping order for some loops
 // !!!TODO: put this in an if statement to only run if the first year of ART has passed
 internal::run_child_art_initiation(time_step, pars, state_curr, state_next, intermediate);
 internal::run_child_art_mortality(time_step, pars, state_curr, state_next, intermediate);
 internal::run_wlhiv_births(time_step, pars, state_curr, state_next, intermediate);


}



} // namespace leapfrog

