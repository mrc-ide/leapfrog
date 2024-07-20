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

template<typename ModelVariant>
constexpr int hBF_coarse = StateSpace<ModelVariant>().children.hBF_coarse;

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
  TensorFixedSize <real_type, Sizes<hAG<ModelVariant>>> hivpop_ha;


  // ART initiation and allocation calculations
  real_type Xart_15plus;

  TensorFixedSize <real_type, Sizes<hDS<ModelVariant>, hAG<ModelVariant>>> artelig_hahm;
  TensorFixedSize <real_type, Sizes<hDS<ModelVariant>>> artelig_hm;
  real_type Xartelig_15plus;

  TensorFixedSize <real_type, Sizes<hDS<ModelVariant>>> expect_mort_artelig_hm;
  real_type expect_mort_artelig15plus;

  TensorFixedSize <real_type, Sizes<hDS<ModelVariant>>> artinit_hm;
  

  real_type cd4mx_scale;
  real_type artpop_hahm;
  real_type deaths;
  int everARTelig_idx;
  int cd4elig_idx;
  real_type p_infections_a;
  real_type p_infections_ha;
  real_type deaths_art;

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
    hivpop_ha.setZero();
    cd4mx_scale = 1.0;
    deaths = 0.0;
    everARTelig_idx = 0;
    cd4elig_idx = 0;
    p_infections_a = 0.0;
    p_infections_ha = 0.0;
    deaths_art = 0.0;
    //
    // ART initiation and allocation
    artelig_hahm.setZero();
    artelig_hm.setZero();
    Xart_15plus = 0.0;
    Xartelig_15plus = 0.0;
    expect_mort_artelig_hm.setZero();
    expect_mort_artelig15plus = 0.0;
    artinit_hm.setZero();
    //
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
  TensorFixedSize <real_type, Sizes<hTS<ChildModel>, hDS<ChildModel>, NS<ChildModel>>> age15_art_pop;
  TensorFixedSize <real_type, Sizes<hDS<ChildModel>, hcTT<ChildModel>, hAG<ChildModel>, NS<ChildModel>>> hc_posthivmort;
  TensorFixedSize <real_type, Sizes<hDS<ChildModel>, hcTT<ChildModel>, hAG<ChildModel>, NS<ChildModel>>> hc_grad;
  TensorFixedSize <real_type, Sizes<hDS<ChildModel>, hcTT<ChildModel>, hAG<ChildModel>, NS<ChildModel>>> hc_ctx_need;
  TensorFixedSize <real_type, Sizes<hDS<ChildModel>, hAG<ChildModel>, NS<ChildModel>>> eligible;
  TensorFixedSize <real_type, Sizes<4>> unmet_need;
  TensorFixedSize <real_type, Sizes<4>> total_need;
  TensorFixedSize <real_type, Sizes<4>> on_art;
  real_type retained;
  TensorFixedSize <real_type, Sizes<4>> total_art_last_year;
  TensorFixedSize <real_type, Sizes<4>> total_art_this_year;
  real_type hc_death_rate;
  TensorFixedSize <real_type, Sizes<hDS<ChildModel>, hcTT<ChildModel>, hAG<ChildModel>, NS<ChildModel>>> hc_art_grad;
  TensorFixedSize <real_type, Sizes<4>> hc_art_scalar;
  TensorFixedSize <real_type, Sizes<4>> hc_initByAge;
  TensorFixedSize <real_type, Sizes<4>> hc_adj;
  TensorFixedSize <real_type, Sizes<4>> hc_art_deaths;
  real_type asfr_sum;
  real_type births_sum;
  real_type nHIVcurr;
  real_type nHIVlast;
  real_type df;
  real_type prev;
  real_type birthsCurrAge;
  real_type birthsHE;
  real_type births_HE_15_24;
  //started here
  TensorFixedSize <real_type, Sizes<hAG<ChildModel>, NS<ChildModel>>> p_hiv_neg_pop;
  real_type sumARV;
  real_type need_PMTCT;
  TensorFixedSize <real_type, Sizes<hPS<ChildModel>>> PMTCT_coverage;
  real_type OnPMTCT;
  real_type num_wlhiv_lt200;
  real_type num_wlhiv_200to350;
  real_type num_wlhiv_gte350;
  real_type num_wlhiv;
  real_type prop_wlhiv_lt200;
  real_type prop_wlhiv_200to350;
  real_type prop_wlhiv_gte350;
  real_type prop_wlhiv_lt350;
  real_type excessratio;
  real_type excessratio_bf;
  real_type optA_transmission_rate;
  real_type optB_transmission_rate;
  real_type optA_bf_transmission_rate;
  real_type optB_bf_transmission_rate;
  real_type retained_on_ART;
  real_type retained_started_ART;
  real_type perinatal_transmission_rate;
  real_type receiving_PMTCT;
  real_type no_PMTCT;
  real_type perinatal_transmission_rate_bf_calc;
  real_type age_weighted_hivneg;
  real_type age_weighted_infections;
  real_type incidence_rate_wlhiv;
  real_type perinatal_transmission_from_incidence;
  real_type bf_at_risk;
  real_type bf_incident_hiv_transmission_rate;
  real_type percent_no_treatment;
  real_type percent_on_treatment;
  TensorFixedSize <real_type, Sizes<hBF_coarse<ChildModel>>> bf_transmission_rate;
  real_type ctx_coverage;
  real_type need_cotrim;


  ChildModelIntermediateData() {};

  void reset() {
    age15_hiv_pop.setZero();
    age15_art_pop.setZero();
    hc_posthivmort.setZero();
    hc_grad.setZero();
    hc_ctx_need.setZero();
    eligible.setZero();
    unmet_need.setZero();
    total_need.setZero();
    on_art.setZero();
    retained = 0.0;
    total_art_last_year.setZero();
    total_art_this_year.setZero();
    hc_death_rate = 0.0;
    hc_art_grad.setZero();
    hc_art_scalar.setZero();
    hc_initByAge.setZero();
    hc_adj.setZero();
    hc_art_deaths.setZero();
    asfr_sum = 0.0;
    births_sum = 0.0;
    nHIVcurr = 0.0;
    nHIVlast = 0.0;
    df = 0.0;
    prev = 0.0;
    birthsCurrAge = 0.0;
    birthsHE = 0.0;
    births_HE_15_24 = 0.0;
    p_hiv_neg_pop.setZero();
    sumARV = 0.0;
    need_PMTCT = 0.0;
    PMTCT_coverage.setZero();
    OnPMTCT = 0.0;
    num_wlhiv_lt200 = 0.0;
    num_wlhiv_200to350 = 0.0;
    num_wlhiv_gte350 = 0.0;
    num_wlhiv = 0.0;
    prop_wlhiv_lt200 = 0.0;
    prop_wlhiv_200to350 = 0.0;
    prop_wlhiv_gte350 = 0.0;
    prop_wlhiv_lt350 = 0.0;
    excessratio = 0.0;
    excessratio_bf = 0.0;
    optA_transmission_rate = 0.0;
    optB_transmission_rate = 0.0;
    optA_bf_transmission_rate = 0.0;
    optB_bf_transmission_rate = 0.0;
    retained_on_ART = 0.0;
    retained_started_ART = 0.0;
    perinatal_transmission_rate = 0.0;
    receiving_PMTCT = 0.0;
    no_PMTCT = 0.0;
    perinatal_transmission_rate_bf_calc = 0.0;
    age_weighted_hivneg = 0.0;
    age_weighted_infections = 0.0;
    incidence_rate_wlhiv = 0.0;
    perinatal_transmission_from_incidence = 0.0;
    bf_at_risk = 0.0;
    bf_incident_hiv_transmission_rate = 0.0;
    percent_no_treatment = 0.0;
    percent_on_treatment = 0.0;
    bf_transmission_rate.setZero();
    ctx_coverage = 0.0;
    need_cotrim = 0.0;
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
