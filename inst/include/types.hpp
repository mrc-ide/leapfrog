#pragma once

#include <unsupported/Eigen/CXX11/Tensor>

namespace leapfrog {

  const int num_genders = 2;
  const int age_groups_pop = 81;
  
  const int pIDX_HIVADULT = 15;
  const int pIDX_FERT = 15;
  const int pAG_FERT = 35;
  // const int hAG_COARSE = 9;
  // const int hAG_FULL = 66;
  const int age_groups_hiv = 66;
  const int age_groups_hiv_15plus = 66;
  const int disease_stages = 7;
  const int treatment_stages = 3;  

  
  // JE ADDED
  using Eigen::Sizes;
  using Eigen::TensorFixedSize;
  
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

  // template<typename real_type>
  // using Tensor2 = Eigen::Tensor<real_type, 2>;

  // template<typename real_type>
  // using Tensor3 = Eigen::Tensor<real_type, 3>;

  // template<typename real_type>
  // using Tensor4 = Eigen::Tensor<real_type, 4>;

template<typename real_type>
struct Parameters {
  int num_genders;
  // Default 81 for ages 0 to 80+
  int age_groups_pop;
  // First index of population eligible for fertility
  int fertility_first_age_group;
  // Number of ages eligible for fertility
  int age_groups_fert;
  // Numer of age groups in HIV population
  int age_groups_hiv;
  // Numer of age groups in HIV 15+ population
  int age_groups_hiv_15plus;
  // Number of HIV disease stages
  int disease_stages;
  // First index of HIV population to model as adult
  int hiv_adult_first_age_group;
  // Number of HIV ART treatment stages
  int treatment_stages;
  // Time step to start ART treatment
  int time_art_start;
  // Index of the youngest age that is reflected in the adult incidence input
  int adult_incidence_first_age_group;

  int pAG_INCIDPOP;
  // Number of time steps per year in the HIV projection
  int hiv_steps_per_year;
  // Difference in time for each hiv time step in HIV projection
  double dt;

  int scale_cd4_mortality;

  int hIDX_15PLUS;

  real_type art_alloc_mxweight;

  // Number of years in each HIV age group
  TensorMap1<int> age_groups_hiv_span;

  // Incidence rate at each time step
  TensorMap1<real_type> incidence_rate;

  TensorMap2<real_type> base_pop;
  TensorMap3<real_type> survival;
  TensorMap3<real_type> net_migration;
  TensorMap2<real_type> age_sex_fertility_ratio;
  TensorMap2<real_type> births_sex_prop;
  TensorMap3<real_type> incidence_relative_risk_age;
  TensorMap1<real_type> incidence_relative_risk_sex;
  TensorMap3<real_type> cd4_mortality;
  TensorMap3<real_type> cd4_progression;
  Tensor1<int> artcd4elig_idx;
  TensorMap3<real_type> cd4_initdist;
  TensorMap1<int> hiv_age_groups_span;
  TensorMap4<real_type> art_mortality;
  TensorMap2<real_type> artmx_timerr;
  Tensor1<real_type> h_art_stage_dur;
  TensorMap1<real_type> art_dropout;
  TensorMap2<real_type> art15plus_num;
  TensorMap2<int> art15plus_isperc;
};

template<typename real_type>
struct State {
  // Tensor2<real_type> total_population;
  // Tensor2<real_type> natural_deaths;
  // Tensor2<real_type> hiv_population;
  // Tensor2<real_type> hiv_natural_deaths;
  // Tensor3<real_type> hiv_strat_adult;
  // Tensor4<real_type> art_strat_adult;
  // real_type births;
  // Tensor3<real_type> aids_deaths_no_art;
  // Tensor2<real_type> infections;
  // Tensor4<real_type> aids_deaths_art;
  // Tensor3<real_type> art_initiation;
  // Tensor2<real_type> hiv_deaths;

  TensorFixedSize<real_type, Sizes<age_groups_pop, num_genders>> total_population;
  TensorFixedSize<real_type, Sizes<age_groups_pop, num_genders>> natural_deaths;
  TensorFixedSize<real_type, Sizes<age_groups_pop, num_genders>> hiv_population;
  TensorFixedSize<real_type, Sizes<age_groups_pop, num_genders>> hiv_natural_deaths;
  TensorFixedSize<real_type, Sizes<disease_stages, age_groups_hiv, num_genders>> hiv_strat_adult;
  TensorFixedSize<real_type, Sizes<treatment_stages,
				   disease_stages,
				   age_groups_hiv,
				   num_genders>> art_strat_adult;
  real_type births;
  TensorFixedSize<real_type, Sizes<disease_stages, age_groups_hiv, num_genders>> aids_deaths_no_art;
  TensorFixedSize<real_type, Sizes<age_groups_pop, num_genders>> infections;
  TensorFixedSize<real_type, Sizes<treatment_stages, disease_stages, age_groups_hiv, num_genders>> aids_deaths_art;
  TensorFixedSize<real_type, Sizes<disease_stages, age_groups_hiv, num_genders>> art_initiation;
  TensorFixedSize<real_type, Sizes<age_groups_pop, num_genders>> hiv_deaths;

  State(int age_groups_pop,
        int num_genders,
        int disease_stages,
        int age_groups_hiv,
        int treatment_stages)
      // : total_population(age_groups_pop, num_genders),
      //   natural_deaths(age_groups_pop, num_genders),
      //   hiv_population(age_groups_pop, num_genders),
      //   hiv_natural_deaths(age_groups_pop, num_genders),
      //   hiv_strat_adult(disease_stages, age_groups_hiv, num_genders),
      //   art_strat_adult(treatment_stages,
      //                   disease_stages,
      //                   age_groups_hiv,
      //                   num_genders),
      //   aids_deaths_no_art(disease_stages, age_groups_hiv, num_genders),
      //   infections(age_groups_pop, num_genders),
      //   aids_deaths_art(treatment_stages, disease_stages, age_groups_hiv, num_genders),
      //   art_initiation(disease_stages, age_groups_hiv, num_genders),
      //   hiv_deaths(age_groups_pop, num_genders)
      : total_population(),
        natural_deaths(),
        hiv_population(),
        hiv_natural_deaths(),
        hiv_strat_adult(),
        art_strat_adult(),
        aids_deaths_no_art(),
        infections(),
        aids_deaths_art(),
        art_initiation(),
        hiv_deaths()
  {}

  void reset() {
    total_population.setZero();
    natural_deaths.setZero();
    hiv_population.setZero();
    hiv_natural_deaths.setZero();
    hiv_strat_adult.setZero();
    art_strat_adult.setZero();
    aids_deaths_no_art.setZero();
    infections.setZero();
    aids_deaths_art.setZero();
    art_initiation.setZero();
    hiv_deaths.setZero();
    births = 0;
  }
};

namespace internal {

  const int NG = 2;
  const int num_genders = 2;
  const int pAG = 81;
  const int age_groups_pop = 81;
  
  const int pIDX_HIVADULT = 15;
  const int pIDX_FERT = 15;
  const int pAG_FERT = 35;
  // const int hAG_COARSE = 9;
  // const int hAG_FULL = 66;
  const int hAG = 66;
  const int age_groups_hiv = 66;
  const int age_groups_hiv_15plus = 66;
  const int hDS = 7;
  const int disease_stages = 7;
  const int hTS = 3;
  
const int MALE = 0;
const int FEMALE = 1;
const int ART0MOS = 0;

template<typename real_type>
struct IntermediateData {
  // Tensor2<real_type> migration_rate;
  TensorFixedSize<real_type, Sizes<age_groups_pop, num_genders>> migration_rate;
  // Tensor2<real_type> hiv_net_migration;
  TensorFixedSize<real_type, Sizes<age_groups_pop, num_genders>> hiv_net_migration;
  // Tensor2<real_type> hiv_population_coarse_ages;
  TensorFixedSize<real_type, Sizes<age_groups_hiv, num_genders>> hiv_population_coarse_ages;
  // Tensor2<real_type> hiv_age_up_prob;
  TensorFixedSize<real_type, Sizes<age_groups_hiv, num_genders>> hiv_age_up_prob;
  // Tensor2<real_type> hiv_negative_pop;
  TensorFixedSize<real_type, Sizes<age_groups_pop, num_genders>> hiv_negative_pop;
  // Tensor2<real_type> infections_ts;
  TensorFixedSize<real_type, Sizes<age_groups_pop, num_genders>> infections_ts;
  // Tensor1<real_type> incidence_rate_sex;
  TensorFixedSize<real_type, Sizes<num_genders>> incidence_rate_sex;
  // Tensor1<real_type> hiv_neg_aggregate;
  TensorFixedSize<real_type, Sizes<num_genders>> hiv_neg_aggregate;
  // Tensor1<real_type> Xhivn_incagerr;
  TensorFixedSize<real_type, Sizes<num_genders>> Xhivn_incagerr;
  // Tensor2<real_type> hiv_deaths_age_sex;
  TensorFixedSize<real_type, Sizes<age_groups_hiv, num_genders>> hiv_deaths_age_sex;
  // Tensor3<real_type> grad;
  TensorFixedSize<real_type, Sizes<hDS, hAG, NG>> grad;
  // Tensor4<real_type> gradART;
  TensorFixedSize<real_type, Sizes<hTS, hDS, hAG, NG>> gradART;  
  // Tensor2<real_type> artelig_hahm;
  TensorFixedSize<real_type, Sizes<disease_stages, age_groups_hiv_15plus>> artelig_hahm;
  // Tensor1<real_type> hivpop_ha;
  TensorFixedSize<real_type, Sizes<age_groups_hiv>> hivpop_ha;
  real_type cd4mx_scale;
  real_type artpop_hahm;
  real_type deaths;
  int everARTelig_idx;
  int cd4elig_idx;
  real_type infections_a;
  real_type infections_ha;
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

  IntermediateData(int age_groups_pop,
		   int age_groups_hiv,
		   int num_genders,
		   int disease_stages,
                   int treatment_stages,
                   int age_groups_hiv_15plus)
    : // migration_rate(age_groups_pop, num_genders),
    migration_rate(),
    // hiv_net_migration(age_groups_pop, num_genders),
    //     hiv_population_coarse_ages(age_groups_hiv, num_genders),
    //     hiv_age_up_prob(age_groups_hiv, num_genders),
    //     hiv_negative_pop(age_groups_pop, num_genders),
    hiv_net_migration(),
        hiv_population_coarse_ages(),
        hiv_age_up_prob(),
        hiv_negative_pop(),
        // infections_ts(age_groups_pop, num_genders),
        // incidence_rate_sex(num_genders),
        // hiv_neg_aggregate(num_genders),
        // Xhivn_incagerr(num_genders),
        // hiv_deaths_age_sex(age_groups_hiv, num_genders),
        infections_ts(),
        incidence_rate_sex(),
        hiv_neg_aggregate(),
        Xhivn_incagerr(),
        hiv_deaths_age_sex(),    
        // grad(disease_stages, age_groups_hiv, num_genders),
        // gradART(treatment_stages, disease_stages, age_groups_hiv, num_genders),
	grad(),
	gradART(),	
        // artelig_hahm(disease_stages, age_groups_hiv_15plus),
        // hivpop_ha(age_groups_hiv),
        artelig_hahm(),
        hivpop_ha(),
        cd4mx_scale(1.0),
        artpop_hahm(0.0),
        deaths(0.0),
        everARTelig_idx(0),
        cd4elig_idx(0),
        infections_a(0.0),
        infections_ha(0.0),
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
    hiv_population_coarse_ages.setZero();
    hiv_age_up_prob.setZero();
    hiv_negative_pop.setZero();
    infections_ts.setZero();
    hiv_neg_aggregate.setZero();
    Xhivn_incagerr.setZero();
    incidence_rate_sex.setZero();
    hiv_deaths_age_sex.setZero();
    grad.setZero();
    gradART.setZero();
    artelig_hahm.setZero();
    hivpop_ha.setZero();
    cd4mx_scale = 1.0;
    deaths = 0.0;
    everARTelig_idx = 0;
    cd4elig_idx = 0;
    infections_a = 0.0;
    infections_ha = 0.0;
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
