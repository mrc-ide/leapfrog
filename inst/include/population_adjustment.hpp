#pragma once

#include <iostream>
#include "types.hpp"

namespace leapfrog {


template<typename ModelVariant, typename real_type>
void calc_hiv_negative_population(
    int time_step,
    const Parameters<ModelVariant, real_type> &pars,
    const State<ModelVariant, real_type> &state_curr,
    State<ModelVariant, real_type> &state_next,
    internal::IntermediateData<ModelVariant, real_type> &intermediate) {
  constexpr auto ss = StateSpace<ModelVariant>().base;

  for (int g = 0; g < ss.NS; ++g) {
    for (int a = 0; a < ss.pAG; ++a) {
        intermediate.base.hivneg_pop(a, g) = state_curr.base.p_total_pop(a, g) - state_curr.base.p_hiv_pop(a, g);
    }
  }

}

template<typename ModelVariant, typename real_type>
void run_base_population_adjustment(
    int time_step,
    const Parameters<ModelVariant, real_type> &pars,
    const State<ModelVariant, real_type> &state_curr,
    State<ModelVariant, real_type> &state_next,
    internal::IntermediateData<ModelVariant, real_type> &intermediate) {
  constexpr auto ss = StateSpace<ModelVariant>().base;
  const auto demog = pars.base.demography;
  // internal::calc_hiv_negative_population(time_step, pars, state_curr, state_next, intermediate);
  for (int g = 0; g < ss.NS; ++g) {
    for (int a = 0; a < ss.pAG; ++a) {
      intermediate.base.hivneg_pop(a, g) = state_curr.base.p_total_pop(a, g) - state_curr.base.p_hiv_pop(a, g);
    }
  }
  for(int g = 0; g < ss.NS; ++g){
    int a = 0;
    for(int ha = 0; ha < ss.hAG; ha++){
      for(int i = 0; i < ss.hAG_span[ha]; i++){
        //base_pop <== targetpop
        //hivpop <== pop?
        intermediate.base.hivadj_ha += state_next.base.p_hiv_pop(a,g);

        intermediate.base.popadjrate_a = demog.base_pop(a,g,time_step) / (intermediate.base.hivneg_pop(a, g) + state_next.base.p_hiv_pop(a, g));

        intermediate.base.hpopadj_a = (intermediate.base.popadjrate_a - 1.0) * state_curr.base.p_hiv_pop(a, g);
        intermediate.base.popadj_ha += intermediate.base.hpopadj_a;
     //   state_curr.base.p_hiv_pop(a, g) += intermediate.base.hpopadj_a;
        a++;
      }

    // population adjustment for hivpop
        if(intermediate.base.hivadj_ha > 0 ){
        intermediate.base.popadjrate_ha = intermediate.base.popadj_ha / intermediate.base.hivadj_ha;
        }else{
         intermediate.base.popadjrate_ha = 0.0;
        }

  if(a >= pars.base.options.p_idx_fertility_first){
    for(int hm = 0; hm < ss.hDS; ++hm){
      state_next.base.h_hiv_adult(hm, ha, g) *= 1+intermediate.base.popadjrate_ha;
      if(time_step >= pars.base.options.ts_art_start)
        for(int hu = 0; hu < ss.hTS; ++hu)
          state_next.base.h_art_adult(hu, hm, ha, g) *= 1+intermediate.base.popadjrate_ha;
    } // loop over hm
  }

    } // loop over ha
  } // loop over g

}

template<typename ModelVariant, typename real_type>
void run_child_population_adjustment(
    int time_step,
    const Parameters<ModelVariant, real_type> &pars,
    const State<ModelVariant, real_type> &state_curr,
    State<ModelVariant, real_type> &state_next,
    internal::IntermediateData<ModelVariant, real_type> &intermediate) {
  constexpr auto ss = StateSpace<ModelVariant>().base;
  constexpr auto hc_ss = StateSpace<ModelVariant>().children;
  const auto demog = pars.base.demography;
  // internal::calc_hiv_negative_population(time_step, pars, state_curr, state_next, intermediate);
  for (int g = 0; g < ss.NS; ++g) {
    for (int a = 0; a < ss.pAG; ++a) {
      intermediate.base.hivneg_pop(a, g) = state_curr.base.p_total_pop(a, g) - state_curr.base.p_hiv_pop(a, g);
    }
  }

  for(int g = 0; g < ss.NS; ++g){
    for(int a = 0; a < pars.base.options.p_idx_fertility_first; a++){

        intermediate.children.hivadj_ha += state_next.base.p_hiv_pop(a,g);
        intermediate.children.popadjrate_a = demog.base_pop(a,g,time_step) / (intermediate.base.hivneg_pop(a, g) + state_next.base.p_hiv_pop(a, g));
        intermediate.children.hpopadj_a = (intermediate.children.popadjrate_a - 1.0) * state_curr.base.p_hiv_pop(a, g);
        intermediate.children.popadj_ha += intermediate.children.hpopadj_a;
   //     state_curr.base.p_hiv_pop(a, g) += intermediate.children.hpopadj_a;

      // population adjustment for hivpop
      if(intermediate.children.hivadj_ha > 0 ){
        intermediate.children.popadjrate_ha = intermediate.children.popadj_ha / intermediate.children.hivadj_ha;
      }else{
        intermediate.children.popadjrate_ha = 0.0;
      }

      if(a < hc_ss.hc1AG){
        for(int hm = 0; hm < hc_ss.hc1DS; ++hm){
          for (int cat = 0; cat < hc_ss.hcTT; ++cat) {
          state_next.children.hc1_hiv_pop(hm, cat, a, g) *= 1+intermediate.children.popadjrate_ha;
          }
          if(time_step >= pars.base.options.ts_art_start)
            for(int hu = 0; hu < ss.hTS; ++hu){
              state_next.children.hc1_art_pop(hu, hm, a, g) *= 1+intermediate.children.popadjrate_ha;
            }
        } // loop over hm
      }else{
        for(int hm = 0; hm < hc_ss.hc2DS; ++hm){
          for (int cat = 0; cat < hc_ss.hcTT; ++cat) {
            state_next.children.hc2_hiv_pop(hm, cat, a, g) *= 1+intermediate.children.popadjrate_ha;
          }
          if(time_step >= pars.base.options.ts_art_start)
            for(int hu = 0; hu < ss.hTS; ++hu){
              state_next.children.hc2_art_pop(hu, hm, a, g) *= 1+intermediate.children.popadjrate_ha;
            }
        } // loop over hm
      }
    } // loop over ha
  } // loop over g

}

}
