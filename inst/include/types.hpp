#pragma once

#include <unsupported/Eigen/CXX11/Tensor>
#include "state_space.hpp"

namespace leapfrog {

template<typename real_type>
using TensorMap1 = Eigen::TensorMap <Eigen::Tensor<real_type, 1>>;

template<typename real_type>
using TensorMap2 = Eigen::TensorMap <Eigen::Tensor<real_type, 2>>;

template<typename real_type>
using TensorMap3 = Eigen::TensorMap <Eigen::Tensor<real_type, 3>>;

template<typename real_type>
using TensorMap4 = Eigen::TensorMap <Eigen::Tensor<real_type, 4>>;

template<typename real_type>
using Tensor1 = Eigen::Tensor<real_type, 1>;

template<typename real_type>
using Tensor2 = Eigen::Tensor<real_type, 2>;

template<typename real_type>
using Tensor3 = Eigen::Tensor<real_type, 3>;

template<typename real_type>
using Tensor4 = Eigen::Tensor<real_type, 4>;

template<typename real_type>
using Tensor5 = Eigen::Tensor<real_type, 5>;

template<HivAgeStratification S>
constexpr int NS = StateSpace<S>::NS;

template<HivAgeStratification S>
constexpr int pAG = StateSpace<S>::pAG;

template<HivAgeStratification S>
constexpr int hAG = StateSpace<S>::hAG;

template<HivAgeStratification S>
constexpr int hDS = StateSpace<S>::hDS;

template<HivAgeStratification S>
constexpr int hTS = StateSpace<S>::hTS;

template<HivAgeStratification S>
constexpr int hC2_disease_stages = StateSpace<S>::hC2_disease_stages;

template<HivAgeStratification S>
constexpr int hTM = StateSpace<S>::hTM;

template<HivAgeStratification S>
constexpr int hPS = StateSpace<S>::hPS;

template<HivAgeStratification S>
constexpr int hBF = StateSpace<S>::hBF;

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
  const int scale_cd4_mortality;
  const int hIDX_15PLUS;
  const real_type initiation_mortality_weight;
  const bool run_child_model;

  Options(int hts_per_year,
          int ts_art_start,
          int hAG_15plus,
          int scale_cd4_mortality,
          real_type initiation_mortality_weight,
          bool run_child_model) :
      hts_per_year(hts_per_year),
      dt(1.0 / hts_per_year),
      p_idx_fertility_first(15),
      p_fertility_age_groups(35),
      p_idx_hiv_first_adult(15),
      adult_incidence_first_age_group(15),
      pAG_INCIDPOP(35),
      ts_art_start(ts_art_start),
      hAG_15plus(hAG_15plus),
      scale_cd4_mortality(scale_cd4_mortality),
      hIDX_15PLUS(0),
      initiation_mortality_weight(initiation_mortality_weight),
      run_child_model(run_child_model) {}
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
};

template<typename real_type>
struct Art {
  Tensor1<int> idx_hm_elig;
  TensorMap4<real_type> mortality;
  TensorMap2<real_type> mortaility_time_rate_ratio;
  Tensor1<real_type> h_art_stage_dur;
  TensorMap1<real_type> dropout;
  TensorMap2<real_type> adults_on_art;
  TensorMap2<int> adults_on_art_is_percent;
};

template<typename real_type>
struct Children {
  TensorMap1<real_type> hc_nosocomial;
  TensorMap1<real_type> hc1_cd4_dist;
  TensorMap2<real_type> hc_cd4_transition;
};

template<typename real_type>
struct Parameters {
  Options<real_type> options;
  Demography<real_type> demography;
  Incidence<real_type> incidence;
  NaturalHistory<real_type> natural_history;
  Art<real_type> art;
  Children<real_type> children;
};

namespace {
using Eigen::Sizes;
using Eigen::TensorFixedSize;
}

template<HivAgeStratification S, typename real_type>
struct State {
  TensorFixedSize <real_type, Sizes<pAG<S>, NS<S>>> p_total_pop;
  TensorFixedSize <real_type, Sizes<pAG<S>, NS<S>>> p_total_pop_natural_deaths;
  TensorFixedSize <real_type, Sizes<pAG<S>, NS<S>>> p_hiv_pop;
  TensorFixedSize <real_type, Sizes<pAG<S>, NS<S>>> p_hiv_pop_natural_deaths;
  TensorFixedSize <real_type, Sizes<hDS<S>, hAG<S>, NS<S>>> h_hiv_adult;
  TensorFixedSize <real_type, Sizes<hTS<S>, hDS<S>, hAG<S>, NS<S>>>
      h_art_adult;
  real_type births;
  TensorFixedSize <real_type, Sizes<hDS<S>, hAG<S>, NS<S>>> h_hiv_deaths_no_art;
  TensorFixedSize <real_type, Sizes<pAG<S>, NS<S>>> p_infections;
  TensorFixedSize <real_type, Sizes<hTS<S>, hDS<S>, hAG<S>, NS<S>>>
    h_hiv_deaths_art;
  TensorFixedSize <real_type, Sizes<hDS<S>, hAG<S>, NS<S>>> h_art_initiation;
  TensorFixedSize <real_type, Sizes<pAG<S>, NS<S>>> p_hiv_deaths;
  TensorFixedSize <real_type, Sizes<hDS<S>, hTM<S>, pAG<S>, NS<S>>> hc_hiv_pop;

  State() {}

  void reset() {
    p_total_pop.setZero();
    p_total_pop_natural_deaths.setZero();
    p_hiv_pop.setZero();
    p_hiv_pop_natural_deaths.setZero();
    h_hiv_adult.setZero();
    h_art_adult.setZero();
    h_hiv_deaths_no_art.setZero();
    p_infections.setZero();
    h_hiv_deaths_art.setZero();
    h_art_initiation.setZero();
    p_hiv_deaths.setZero();
    hc_hiv_pop.setZero();
    births = 0;
  }
};

namespace internal {

const int MALE = 0;
const int FEMALE = 1;
const int ART0MOS = 0;

template<HivAgeStratification S, typename real_type>
struct IntermediateData {
  TensorFixedSize <real_type, Sizes<pAG<S>, NS<S>>> migration_rate;
  TensorFixedSize <real_type, Sizes<pAG<S>, NS<S>>> hiv_net_migration;
  TensorFixedSize <real_type, Sizes<hAG<S>, NS<S>>> p_hiv_pop_coarse_ages;
  TensorFixedSize <real_type, Sizes<hAG<S>, NS<S>>> hiv_age_up_prob;
  TensorFixedSize <real_type, Sizes<pAG<S>, NS<S>>> hiv_negative_pop;
  TensorFixedSize <real_type, Sizes<pAG<S>, NS<S>>> p_infections_ts;
  TensorFixedSize <real_type, Sizes<NS<S>>> rate_sex;
  TensorFixedSize <real_type, Sizes<NS<S>>> hiv_neg_aggregate;
  TensorFixedSize <real_type, Sizes<NS<S>>> Xhivn_incagerr;
  TensorFixedSize <real_type, Sizes<hAG<S>, NS<S>>> p_hiv_deaths_age_sex;
  TensorFixedSize <real_type, Sizes<hDS<S>, hAG<S>, NS<S>>> grad;
  TensorFixedSize <real_type, Sizes<hTS<S>, hDS<S>, hAG<S>, NS<S>>> gradART;
  Tensor2<real_type> artelig_hahm;
  TensorFixedSize <real_type, Sizes<hAG<S>>> hivpop_ha;
  TensorFixedSize <real_type, Sizes<hDS<S>, NS<S>>> age15_hiv_pop;
  real_type cd4mx_scale;
  real_type artpop_hahm;
  real_type deaths;
  int everARTelig_idx;
  int cd4elig_idx;
  real_type p_infections_a;
  real_type p_infections_ha;
  real_type deaths_art;
  real_type Xart_15plus;
  real_type Xartelig_15plus;
  real_type expect_mort_artelig15plus;
  int anyelig_idx;
  real_type artnum_hts;
  real_type artcov_hts;
  real_type curr_coverage;
  real_type artinit_hts;
  real_type artinit_hahm;
  real_type hivqx_ha;
  real_type hivdeaths_a;

  IntermediateData(int hAG_15plus)
      :
      artelig_hahm(hDS<S>, hAG_15plus),
      cd4mx_scale(1.0),
      artpop_hahm(0.0),
      deaths(0.0),
      everARTelig_idx(0),
      cd4elig_idx(0),
      p_infections_a(0.0),
      p_infections_ha(0.0),
      deaths_art(0.0),
      Xart_15plus(0.0),
      Xartelig_15plus(0.0),
      expect_mort_artelig15plus(0.0),
      anyelig_idx(0),
      artnum_hts(0.0),
      artcov_hts(0.0),
      curr_coverage(0.0),
      artinit_hts(0.0),
      artinit_hahm(0.0),
      hivqx_ha(0.0),
      hivdeaths_a(0.0) {}

  void reset() {
    migration_rate.setZero();
    hiv_net_migration.setZero();
    p_hiv_pop_coarse_ages.setZero();
    hiv_age_up_prob.setZero();
    hiv_negative_pop.setZero();
    p_infections_ts.setZero();
    hiv_neg_aggregate.setZero();
    Xhivn_incagerr.setZero();
    rate_sex.setZero();
    p_hiv_deaths_age_sex.setZero();
    grad.setZero();
    gradART.setZero();
    artelig_hahm.setZero();
    hivpop_ha.setZero();
    age15_hiv_pop.setZero();
    cd4mx_scale = 1.0;
    deaths = 0.0;
    everARTelig_idx = 0;
    cd4elig_idx = 0;
    p_infections_a = 0.0;
    p_infections_ha = 0.0;
    deaths_art = 0.0;
    Xart_15plus = 0.0;
    Xartelig_15plus = 0.0;
    expect_mort_artelig15plus = 0.0;
    anyelig_idx = 0;
    artnum_hts = 0.0;
    artcov_hts = 0.0;
    curr_coverage = 0.0;
    artinit_hts = 0.0;
    artinit_hahm = 0.0;
    hivqx_ha = 0.0;
    hivdeaths_a = 0.0;
  }
};

}
}
