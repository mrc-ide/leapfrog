#pragma once

#include "types.hpp"
#include "state_space.hpp"
#include "model_variants.hpp"

namespace leapfrog {

template<typename real_type>
struct Options {
  int hts_per_year;
  double dt;
  const int p_idx_fertility_first;
  const int p_fertility_age_groups;
  const int p_idx_hiv_first_adult;
  const int adult_incidence_first_age_group;
  const int pAG_INCIDPOP;
  const int ts_art_start;
  const int hAG_15plus;
  const int hIDX_15PLUS;

  Options(int hts_per_year,
          int ts_art_start,
          int hAG_15plus) :
      hts_per_year(hts_per_year),
      dt(1.0 / hts_per_year),
      p_idx_fertility_first(15),
      p_fertility_age_groups(35),
      p_idx_hiv_first_adult(15),
      adult_incidence_first_age_group(15),
      pAG_INCIDPOP(35),
      ts_art_start(ts_art_start),
      hAG_15plus(hAG_15plus),
      hIDX_15PLUS(0) {}
};

template<typename real_type>
struct Demography {
  TensorMap2<real_type> base_pop;
  TensorMap3<real_type> survival_probability;
  TensorMap3<real_type> net_migration;
  TensorMap2<real_type> age_specific_fertility_rate;
  TensorMap2<real_type> births_sex_prop;
};

template<typename real_type>
struct Incidence {
  TensorMap1<real_type> total_rate;
  TensorMap3<real_type> relative_risk_age;
  TensorMap1<real_type> relative_risk_sex;
};

template<typename real_type>
struct NaturalHistory {
  TensorMap3<real_type> cd4_mortality;
  TensorMap3<real_type> cd4_progression;
  TensorMap3<real_type> cd4_initial_distribution;
  int scale_cd4_mortality;
};

template<typename real_type>
struct Art {
  Tensor1<int> idx_hm_elig;
  TensorMap4<real_type> mortality;
  TensorMap2<real_type> mortaility_time_rate_ratio;
  TensorMap1<real_type> dropout;
  TensorMap2<real_type> adults_on_art;
  TensorMap2<int> adults_on_art_is_percent;
  TensorMap1<real_type> h_art_stage_dur;
  real_type initiation_mortality_weight;
};

template<typename real_type>
struct Children {
  TensorMap1<real_type> hc_nosocomial;
  TensorMap1<real_type> hc1_cd4_dist;
  TensorMap2<real_type> hc_cd4_transition;
  TensorMap3<real_type> hc1_cd4_mort;
  TensorMap3<real_type> hc2_cd4_mort;
  TensorMap1<real_type> hc1_cd4_prog;
  TensorMap1<real_type> hc2_cd4_prog;
  real_type ctx_effect;
  TensorMap1<real_type> ctx_val;
  TensorMap1<real_type> hc_art_elig_age;
  Tensor2<real_type> hc_art_elig_cd4;
  TensorMap3<real_type> hc_art_mort_rr;
  TensorMap3<real_type> hc1_art_mort;
  TensorMap3<real_type> hc2_art_mort;
  TensorMap1<int> hc_art_isperc;
  TensorMap1<real_type> hc_art_val;
  TensorMap2<real_type> hc_art_init_dist;
};


template<typename real_type>
struct BaseModelParameters {
  Options<real_type> options;
  Demography<real_type> demography;
  Incidence<real_type> incidence;
  NaturalHistory<real_type> natural_history;
  Art<real_type> art;
};

template<typename ModelVariant, typename real_type>
struct ChildModelParameters {
};

template<typename real_type>
struct ChildModelParameters<ChildModel, real_type> {
  Children<real_type> children;
};

template<typename ModelVariant, typename real_type>
struct Parameters {
  BaseModelParameters<real_type> base;
  ChildModelParameters<ModelVariant, real_type> children;
};

}
