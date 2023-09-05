#pragma once

#include <unsupported/Eigen/CXX11/Tensor>
#include "state_space.hpp"
#include "model_variants.hpp"

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

template<typename ModelVariant>
constexpr int NS = StateSpace<ModelVariant>().base.NS;

template<typename ModelVariant>
constexpr int pAG = StateSpace<ModelVariant>().base.pAG;

template<typename ModelVariant>
constexpr int hAG = StateSpace<ModelVariant>().base.hAG;

template<typename ModelVariant>
constexpr int hDS = StateSpace<ModelVariant>().base.hDS;

template<typename ModelVariant>
constexpr int hTS = StateSpace<ModelVariant>().base.hTS;

template<typename ModelVariant>
constexpr int hc1DS = StateSpace<ModelVariant>().children.hc1DS;

template<typename ModelVariant>
constexpr int hc1_ageend = StateSpace<ModelVariant>().children.hc1_ageend;

template<typename ModelVariant>
constexpr int hc2_agestart = StateSpace<ModelVariant>().children.hc2_agestart;

template<typename ModelVariant>
constexpr int hc1AG = StateSpace<ModelVariant>().children.hc1AG;

template<typename ModelVariant>
constexpr int hc2AG = StateSpace<ModelVariant>().children.hc2AG;

template<typename ModelVariant>
constexpr int hc2DS = StateSpace<ModelVariant>().children.hc2DS;

template<typename ModelVariant>
constexpr int hcTT = StateSpace<ModelVariant>().children.hcTT;

template<typename ModelVariant>
constexpr int hPS = StateSpace<ModelVariant>().children.hPS;

template<typename ModelVariant>
constexpr int hBF = StateSpace<ModelVariant>().children.hBF;

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
  Tensor1<real_type> h_art_stage_dur;
  TensorMap1<real_type> dropout;
  TensorMap2<real_type> adults_on_art;
  TensorMap2<int> adults_on_art_is_percent;
  real_type initiation_mortality_weight;
};

template<typename real_type>
struct BaseModelParameters {
  Options<real_type> options;
  Demography<real_type> demography;
  Incidence<real_type> incidence;
  NaturalHistory<real_type> natural_history;
  Art<real_type> art;
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

namespace {
using Eigen::Sizes;
using Eigen::TensorFixedSize;
}

template<typename ModelVariant, typename real_type>
struct BaseModelState {
  TensorFixedSize <real_type, Sizes<pAG<ModelVariant>, NS<ModelVariant>>> p_total_pop;
  TensorFixedSize <real_type, Sizes<pAG<ModelVariant>, NS<ModelVariant>>> p_total_pop_natural_deaths;
  TensorFixedSize <real_type, Sizes<pAG<ModelVariant>, NS<ModelVariant>>> p_hiv_pop;
  TensorFixedSize <real_type, Sizes<pAG<ModelVariant>, NS<ModelVariant>>> p_hiv_pop_natural_deaths;
  TensorFixedSize <real_type, Sizes<hDS<ModelVariant>, hAG<ModelVariant>, NS<ModelVariant>>> h_hiv_adult;
  TensorFixedSize <real_type, Sizes<hTS<ModelVariant>, hDS<ModelVariant>, hAG<ModelVariant>, NS<ModelVariant>>>
      h_art_adult;
  real_type births;
  TensorFixedSize <real_type, Sizes<hDS<ModelVariant>, hAG<ModelVariant>, NS<ModelVariant>>> h_hiv_deaths_no_art;
  TensorFixedSize <real_type, Sizes<pAG<ModelVariant>, NS<ModelVariant>>> p_infections;
  TensorFixedSize <real_type, Sizes<hTS<ModelVariant>, hDS<ModelVariant>, hAG<ModelVariant>, NS<ModelVariant>>>
      h_hiv_deaths_art;
  TensorFixedSize <real_type, Sizes<hDS<ModelVariant>, hAG<ModelVariant>, NS<ModelVariant>>> h_art_initiation;
  TensorFixedSize <real_type, Sizes<pAG<ModelVariant>, NS<ModelVariant>>> p_hiv_deaths;

  BaseModelState(const Parameters<ModelVariant, real_type> &pars) {
    constexpr auto ss = StateSpace<ModelVariant>().base;
    for (int g = 0; g < ss.NS; ++g) {
      for (int a = 0; a < ss.pAG; ++a) {
        p_total_pop(a, g) = pars.base.demography.base_pop(a, g);
      }
    }
    p_total_pop_natural_deaths.setZero();
    p_hiv_pop.setZero();
    p_hiv_pop_natural_deaths.setZero();
    h_hiv_adult.setZero();
    h_art_adult.setZero();
    births = 0;
    h_hiv_deaths_no_art.setZero();
    p_infections.setZero();
    h_hiv_deaths_art.setZero();
    h_art_initiation.setZero();
    p_hiv_deaths.setZero();
  }

  void set_initial_state(const Parameters<ModelVariant, real_type> &pars) {
    constexpr auto ss = StateSpace<ModelVariant>().base;
    for (int g = 0; g < ss.NS; ++g) {
      for (int a = 0; a < ss.pAG; ++a) {
        p_total_pop(a, g) = pars.base.demography.base_pop(a, g);
      }
    }
    p_total_pop_natural_deaths.setZero();
    p_hiv_pop.setZero();
    p_hiv_pop_natural_deaths.setZero();
    h_hiv_adult.setZero();
    h_art_adult.setZero();
    births = 0;
    h_hiv_deaths_no_art.setZero();
    p_infections.setZero();
    h_hiv_deaths_art.setZero();
    h_art_initiation.setZero();
    p_hiv_deaths.setZero();
  }

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
    births = 0;
  }
};

template<typename ModelVariant, typename real_type>
struct ChildModelState {
  ChildModelState(const Parameters<ModelVariant, real_type> &pars) {}

  void set_initial_state(const Parameters<ModelVariant, real_type> &pars) {}

  void reset() {}
};

template<typename real_type>
struct ChildModelState<ChildModel, real_type> {
  TensorFixedSize <real_type, Sizes<hc1DS<ChildModel>, hcTT<ChildModel>, hc1AG<ChildModel>, NS<ChildModel>>> hc1_hiv_pop;
  TensorFixedSize <real_type, Sizes<hc2DS<ChildModel>, hcTT<ChildModel>, hc2AG<ChildModel>, NS<ChildModel>>> hc2_hiv_pop;
  TensorFixedSize <real_type, Sizes<hTS<ChildModel>, hc1DS<ChildModel>, hc1AG<ChildModel>, NS<ChildModel>>> hc1_art_pop;
  TensorFixedSize <real_type, Sizes<hTS<ChildModel>, hc2DS<ChildModel>, hc2AG<ChildModel>, NS<ChildModel>>> hc2_art_pop;

  ChildModelState(const Parameters<ChildModel, real_type> &pars) {
    reset();
  }

  void set_initial_state(const Parameters<ChildModel, real_type> &pars) {
    reset();
  }

  void reset() {
    hc1_hiv_pop.setZero();
    hc2_hiv_pop.setZero();
    hc1_art_pop.setZero();
    hc2_art_pop.setZero();
  }
};

template<typename ModelVariant, typename real_type>
struct State {
  BaseModelState<ModelVariant, real_type> base;
  ChildModelState<ModelVariant, real_type> children;

  State(const Parameters<ModelVariant, real_type> &pars) :
      base(pars),
      children(pars) {}

  void set_initial_state(const Parameters<ChildModel, real_type> &pars) {
    base.set_initial_state();
    children.set_initial_state();
  }

  void reset() {
    base.reset();
    children.reset();
  }
};

namespace internal {

const int MALE = 0;
const int FEMALE = 1;
const int ART0MOS = 0;

template<typename ModelVariant, typename real_type>
struct BaseModelIntermediateData {
  TensorFixedSize <real_type, Sizes<pAG<ModelVariant>, NS<ModelVariant>>> migration_rate;
  TensorFixedSize <real_type, Sizes<pAG<ModelVariant>, NS<ModelVariant>>> hiv_net_migration;
  TensorFixedSize <real_type, Sizes<hAG<ModelVariant>, NS<ModelVariant>>> p_hiv_pop_coarse_ages;
  TensorFixedSize <real_type, Sizes<hAG<ModelVariant>, NS<ModelVariant>>> hiv_age_up_prob;
  TensorFixedSize <real_type, Sizes<pAG<ModelVariant>>> hiv_negative_pop;
  TensorFixedSize <real_type, Sizes<pAG<ModelVariant>, NS<ModelVariant>>> p_infections_ts;
  TensorFixedSize <real_type, Sizes<NS<ModelVariant>>> rate_sex;
  TensorFixedSize <real_type, Sizes<NS<ModelVariant>>> hiv_neg_aggregate;
  TensorFixedSize <real_type, Sizes<hAG<ModelVariant>, NS<ModelVariant>>> p_hiv_deaths_age_sex;
  TensorFixedSize <real_type, Sizes<hDS<ModelVariant>, hAG<ModelVariant>, NS<ModelVariant>>> grad;
  TensorFixedSize <real_type, Sizes<hTS<ModelVariant>, hDS<ModelVariant>, hAG<ModelVariant>, NS<ModelVariant>>> gradART;
  Tensor2<real_type> artelig_hahm;
  TensorFixedSize <real_type, Sizes<hAG<ModelVariant>>> hivpop_ha;
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
  real_type Xhivn_incagerr;

  BaseModelIntermediateData(int hAG_15plus)
      :
      artelig_hahm(hDS<ModelVariant>, hAG_15plus),
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
      hivdeaths_a(0.0),
      Xhivn_incagerr(0.0) {}

  void reset() {
    migration_rate.setZero();
    hiv_net_migration.setZero();
    p_hiv_pop_coarse_ages.setZero();
    hiv_age_up_prob.setZero();
    hiv_negative_pop.setZero();
    p_infections_ts.setZero();
    hiv_neg_aggregate.setZero();
    rate_sex.setZero();
    p_hiv_deaths_age_sex.setZero();
    grad.setZero();
    gradART.setZero();
    artelig_hahm.setZero();
    hivpop_ha.setZero();
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
    Xhivn_incagerr = 0.0;
  }
};

template<typename ModelVariant, typename real_type>
struct ChildModelIntermediateData {
  ChildModelIntermediateData() {};

  void reset() {};
};

template<typename real_type>
struct ChildModelIntermediateData<ChildModel, real_type> {
  TensorFixedSize <real_type, Sizes<hDS<ChildModel>, NS<ChildModel>>> age15_hiv_pop;
  TensorFixedSize <real_type, Sizes<hDS<ChildModel>, hcTT<ChildModel>, hAG<ChildModel>, NS<ChildModel>>> hc_posthivmort;
  TensorFixedSize <real_type, Sizes<hDS<ChildModel>, hcTT<ChildModel>, hAG<ChildModel>, NS<ChildModel>>> hc_grad;

  ChildModelIntermediateData() {};

  void reset() {
    age15_hiv_pop.setZero();
    hc_posthivmort.setZero();
    hc_grad.setZero();
  };
};

template<typename ModelVariant, typename real_type>
struct IntermediateData {
  BaseModelIntermediateData<ModelVariant, real_type> base;
  ChildModelIntermediateData<ModelVariant, real_type> children;

  IntermediateData(int hAG_15plus)
      :
      base(hAG_15plus),
      children() {}

  void reset() {
    base.reset();
    children.reset();
  }
};

}
}
