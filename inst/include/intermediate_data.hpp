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

namespace {
using Eigen::Sizes;
using Eigen::TensorFixedSize;
}

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
  TensorFixedSize <real_type, Sizes<hDS<ChildModel>, hcTT<ChildModel>, hAG<ChildModel>, NS<ChildModel>>> hc_art_need;
  TensorFixedSize <real_type, Sizes<hDS<ChildModel>, hcTT<ChildModel>, hAG<ChildModel>, NS<ChildModel>>> hc_art_init;
  real_type hc_art_init_total;
  real_type hc_death_rate;
  TensorFixedSize <real_type, Sizes<hDS<ChildModel>, hcTT<ChildModel>, hAG<ChildModel>, NS<ChildModel>>> hc_art_grad;
  real_type hc_art_scalar;
  real_type hc_initByAge;
  real_type hc_adj;

  ChildModelIntermediateData() {};

  void reset() {
    age15_hiv_pop.setZero();
    hc_posthivmort.setZero();
    hc_grad.setZero();
    hc_art_need.setZero();
    hc_art_init.setZero();
    hc_art_init_total = 0.0;
    hc_death_rate = 0.0;
    hc_art_grad.setZero();
    hc_art_scalar = 0.0;
    hc_initByAge = 0.0;
    hc_adj = 0.0;
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
