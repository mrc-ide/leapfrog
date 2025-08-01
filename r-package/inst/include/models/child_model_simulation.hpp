#pragma once

#include "../options.hpp"
#include "../generated/config_mixer.hpp"

namespace leapfrog {
namespace internal {

template<typename Config>
concept ChildModelSimulationEnabled = RunDemographicProjection<Config> && RunHivSimulation<Config> && RunChildModel<Config>;

template<typename Config>
struct ChildModelSimulation {
  ChildModelSimulation(...) {};
};

template<ChildModelSimulationEnabled Config>
struct ChildModelSimulation<Config> {
  using real_type = typename Config::real_type;
  using ModelVariant = typename Config::ModelVariant;
  using SS = Config::SS;
  using Pars = Config::Pars;
  using State = Config::State;
  using Intermediate = Config::Intermediate;
  using Args = Config::Args;

  // private members of this struct
  private:
  // state space
  static constexpr int NS = SS::NS;
  static constexpr int pAG = SS::pAG;
  static constexpr int hDS = SS::hDS;
  static constexpr int hTS = SS::hTS;
  static constexpr int hAG = SS::hAG;
  static constexpr auto hAG_span = SS::hAG_span;
  static constexpr int PROJPERIOD_MIDYEAR = SS::PROJPERIOD_MIDYEAR;
  static constexpr int MALE = SS::MALE;
  static constexpr int FEMALE = SS::FEMALE;
  static constexpr int ART0MOS = SS::ART0MOS;
  static constexpr int hcAG_end = SS::hcAG_end;
  static constexpr int hc2_agestart = SS::hc2_agestart;
  static constexpr int hc1_ageend = SS::hc1_ageend;
  static constexpr int hc1DS = SS::hc1DS;
  static constexpr int hcTT = SS::hcTT;
  static constexpr int hc2DS = SS::hc2DS;
  static constexpr int hPS = SS::hPS;
  static constexpr int hBF = SS::hBF;
  static constexpr int hcAG_coarse = SS::hcAG_coarse;
  static constexpr int p_idx_fertility_first = SS::p_idx_fertility_first;
  static constexpr int p_fertility_age_groups = SS::p_fertility_age_groups;
  static constexpr int p_idx_hiv_first_adult = SS::p_idx_hiv_first_adult;
  static constexpr auto hc_age_coarse = SS::hc_age_coarse;
  static constexpr auto hc_age_coarse_cd4 = SS::hc_age_coarse_cd4;
  static constexpr auto hc_cd4_transition = SS::hc_cd4_transition;

  // function args
  int t;
  const Pars& pars;
  const State& state_curr;
  State& state_next;
  Intermediate& intermediate;
  const Options<real_type>& opts;

  // only exposing the constructor and some methods
  public:
  ChildModelSimulation(Args& args):
    t(args.t),
    pars(args.pars),
    state_curr(args.state_curr),
    state_next(args.state_next),
    intermediate(args.intermediate),
    opts(args.opts)
  {};

  void run_child_model_simulation() {
    const auto& p_hc = pars.hc;

    run_child_ageing();

    if (p_hc.mat_prev_input(t)) {
      run_wlhiv_births_input_mat_prev();
    } else {
      run_wlhiv_births();
    }

    adjust_hiv_births();
    add_infections();
    need_for_cotrim();
    cd4_mortality();
    run_child_hiv_mort();
    add_child_grad();

    if (t >= p_hc.hc_art_start) {
      // First calculate who is eligible for treatment and ART initiates
      eligible_for_treatment();
      art_ltfu();
      calc_art_initiates();

      // Before initiating people on to ART, calculate deaths among CLHIV on ART 6-12 months,
      // the progress them to on ART 6-12 months months and calculate deaths for those on ART >12 months
      on_art_mortality(0, 1);
      progress_time_on_art(0, 1);
      on_art_mortality(2, 2);

      // Initiate CLHIV onto ART
      art_initiation_by_age();

      // Calculate on ART mortality among those on ART < 6 months
      on_art_mortality(0, 1);
      // Progress 6 to 12 mo to 12 plus months
      progress_time_on_art(1, 2);

      // Remove lost to follow ups
      apply_ltfu_to_hivpop();
      apply_ltfu_to_artpop();
    }

    nosocomial_infections();
    fill_total_pop_outputs();
  };

  // private methods that we don't want people to call
  private:
  void run_child_ageing() {
    const auto& p_dp = pars.dp;
    const auto& c_hc = state_curr.hc;
    auto& n_hc = state_next.hc;

    for (int s = 0; s < NS; ++s) {
      // less than 5 because there is a cd4 transition between ages 4 and 5
      for (int a = 1; a < hc2_agestart; ++a) {
        for (int hd = 0; hd < hc1DS; ++hd) {
          for (int cat = 0; cat < hcTT; ++cat) {
            n_hc.hc1_hiv_pop(hd, cat, a, s) += c_hc.hc1_hiv_pop(hd, cat, a - 1, s) * p_dp.survival_probability(a, s, t);
          }
          for (int dur = 0; dur < hTS; ++dur) {
            n_hc.hc1_art_pop(dur, hd, a, s) += c_hc.hc1_art_pop(dur, hd, a - 1, s) * p_dp.survival_probability(a, s, t);
          }
        }
      }
    }

    for (int s = 0; s < NS; ++s) {
      for (int hd = 0; hd < hc1DS; ++hd) {
        for (int hd_alt = 0; hd_alt < hc2DS; ++hd_alt) {
          for (int cat = 0; cat < hcTT; ++cat) {
            n_hc.hc2_hiv_pop(hd_alt, cat, 0, s) += c_hc.hc1_hiv_pop(hd, cat, hc1_ageend, s) *
                                                   p_dp.survival_probability(hc2_agestart, s, t) *
                                                   hc_cd4_transition[hd_alt][hd];
          }
          for (int dur = 0; dur < hTS; ++dur) {
            n_hc.hc2_art_pop(dur, hd_alt, 0, s) += c_hc.hc1_art_pop(dur, hd, hc1_ageend, s) *
                                                   p_dp.survival_probability(hc2_agestart, s, t) *
                                                   hc_cd4_transition[hd_alt][hd];
          }
        }
      }
    }

    for (int s = 0; s < NS; ++s) {
      for (int a = (hc2_agestart + 1); a < p_idx_hiv_first_adult; ++a) {
        for (int hd = 0; hd < hc2DS; ++hd) {
          for (int cat = 0; cat < hcTT; ++cat) {
            n_hc.hc2_hiv_pop(hd, cat, a - hc2_agestart, s) += c_hc.hc2_hiv_pop(hd, cat, a - hc2_agestart - 1, s) *
                                                              p_dp.survival_probability(a, s, t);
          }
          for (int dur = 0; dur < hTS; ++dur) {
            n_hc.hc2_art_pop(dur, hd, a - hc2_agestart, s) += c_hc.hc2_art_pop(dur, hd, a - hc2_agestart - 1, s) *
                                                              p_dp.survival_probability(a, s, t);
          }
        }
      }
    }
  };

  void run_wlhiv_births() {
    const auto& p_dp = pars.dp;
    const auto& p_hc = pars.hc;
    const auto& c_ha = state_curr.ha;
    auto& n_ha = state_next.ha;
    auto& n_dp = state_next.dp;
    auto& n_hc = state_next.hc;
    auto& i_hc = intermediate.hc;

    i_hc.asfr_sum = 0.0;
    for (int a = 0; a < p_fertility_age_groups; ++a) {
      i_hc.asfr_sum += p_dp.age_specific_fertility_rate(a, t);
    } // end a

    for (int a = 0; a < p_fertility_age_groups; ++a) {
      i_hc.nHIVcurr = 0.0;
      i_hc.nHIVlast = 0.0;
      i_hc.df = 0.0;
      // Convert from single years to 5-year age group index
      auto a_fert_idx = std::floor(a/5);

      for (int hd = 0; hd < hDS; ++hd) {
        i_hc.nHIVcurr += n_ha.h_hiv_adult(hd, a, 1);
        i_hc.nHIVlast += c_ha.h_hiv_adult(hd, a, 1);
        for (int ht = 0; ht < hTS; ++ht) {
          i_hc.nHIVcurr += n_ha.h_art_adult(ht, hd, a, 1);
          i_hc.nHIVlast += c_ha.h_art_adult(ht, hd, a, 1);
        } // end hTS
      } // end hDS

      i_hc.prev = i_hc.nHIVcurr / n_dp.p_total_pop(a + 15, 1);

      for (int hd = 0; hd < hDS; ++hd) {
        i_hc.df += p_hc.local_adj_factor * p_hc.fert_mult_by_age(a_fert_idx, t) * p_hc.fert_mult_off_art(hd) *
                   (n_ha.h_hiv_adult(hd, a, 1) + c_ha.h_hiv_adult(hd, a, 1)) / 2;
        // women on ART less than 6 months use the off art fertility multiplier
        i_hc.df += p_hc.local_adj_factor * p_hc.fert_mult_by_age(a_fert_idx, t) * p_hc.fert_mult_off_art(hd) *
                   (n_ha.h_art_adult(0, hd, a, 1) + c_ha.h_art_adult(0, hd, a, 1)) / 2;
        for (int ht = 1; ht < hTS; ++ht) {
          i_hc.df += p_hc.local_adj_factor * p_hc.fert_mult_on_art(a_fert_idx) *
                     (n_ha.h_art_adult(ht, hd, a, 1) + c_ha.h_art_adult(ht, hd, a, 1)) / 2;
        } // end hTS
      } // end hDS


      auto midyear_fertileHIV = (i_hc.nHIVcurr + i_hc.nHIVlast) / 2;
      if (i_hc.nHIVcurr > 0) {
        i_hc.df = i_hc.df / midyear_fertileHIV;
      } else {
        i_hc.df = 1;
      }
      i_hc.birthsCurrAge = midyear_fertileHIV * p_hc.total_fertility_rate(t) *
                           i_hc.df / (i_hc.df * i_hc.prev + 1 - i_hc.prev) *
                           p_dp.age_specific_fertility_rate(a, t) / i_hc.asfr_sum ;

      i_hc.birthsHE += i_hc.birthsCurrAge;
      if (a < 9) {
        i_hc.births_HE_15_24 += i_hc.birthsCurrAge;
      }
    } // end a
    n_hc.hiv_births = i_hc.birthsHE;
  };

  void run_wlhiv_births_input_mat_prev() {
    const auto& p_hc = pars.hc;
    auto& n_hc = state_next.hc;

    n_hc.hiv_births = p_hc.mat_hiv_births(t);
  };

  void calc_hiv_negative_pop() {
    auto& n_ha = state_next.ha;
    auto& n_dp = state_next.dp;
    auto& i_hc = intermediate.hc;

    for (int s = 0; s < NS; ++s) {
      for (int a = 0; a < hAG; ++a) {
        i_hc.p_hiv_neg_pop(a, s) = n_dp.p_total_pop(a, s) - n_ha.p_hiv_pop(a, s);
      }// end a
    }// end s
  };

  void adjust_hiv_births() {
    const auto& p_hc = pars.hc;
    auto& n_hc = state_next.hc;

    if (p_hc.abortion(1, t) == 1) {
      n_hc.hiv_births -= n_hc.hiv_births * p_hc.abortion(0, t);
    } else {
      n_hc.hiv_births -=  p_hc.abortion(0, t);
    }
  };

  void convert_PMTCT_num_to_perc() {
    const auto& p_hc = pars.hc;
    auto& n_hc = state_next.hc;
    auto& i_hc = intermediate.hc;

    i_hc.sumARV = 0.0;
    for (int hp = 0; hp < hPS; ++hp) {
      i_hc.sumARV += p_hc.PMTCT(hp, t);
    }

    i_hc.need_PMTCT = std::max(i_hc.sumARV, n_hc.hiv_births);

    i_hc.on_PMTCT = i_hc.sumARV + p_hc.patients_reallocated(t);
    i_hc.on_PMTCT = std::min(i_hc.on_PMTCT, i_hc.need_PMTCT);

    // replace all instances of coverage input as numbers with percentage covered
    if (p_hc.PMTCT_input_is_percent(t)) {
      for (int hp = 0; hp < hPS; ++hp) {
        i_hc.PMTCT_coverage(hp) = p_hc.PMTCT(hp, t) / 100;
      } // end hPS
      i_hc.sumARV = i_hc.sumARV * i_hc.need_PMTCT;
    } else {
      for (int hp = 0; hp < hPS; ++hp) {
        if (i_hc.sumARV == 0) {
          i_hc.PMTCT_coverage(hp) = 0.0;
        } else {
          i_hc.PMTCT_coverage(hp) = p_hc.PMTCT(hp, t) / i_hc.sumARV *
                                    i_hc.on_PMTCT / i_hc.need_PMTCT;
        }
      } // end hPS
    }

    i_hc.PMTCT_coverage(4) = i_hc.PMTCT_coverage(4) * p_hc.PMTCT_dropout(0, t);
    i_hc.PMTCT_coverage(5) = i_hc.PMTCT_coverage(5) * p_hc.PMTCT_dropout(1, t);
  };

  void convert_PMTCT_pre_bf() {
    const auto& p_hc = pars.hc;
    auto& i_hc = intermediate.hc;

    // TODO: Maggie to confirm why Option A/B alt tr aren't used
    for (int hp = 0; hp < hPS; ++hp) {
      i_hc.PMTCT_coverage(hp) *= 1.0 - p_hc.PMTCT_transmission_rate(4, hp, 0);
    } // end hPS
  };

  void calc_wlhiv_cd4_proportion() {
    const auto& p_hc = pars.hc;
    auto& n_ha = state_next.ha;
    auto& i_hc = intermediate.hc;

    // Option A and B were only authorized for women with greater than 350 CD4, so if the percentage of women
    // on option A/B > the proportion of women in this cd4 category, we assume that some must have a cd4 less than 350
    // option AB will be less effective for these women so we adjust for that

    if (p_hc.mat_prev_input(t)) {
      i_hc.prop_wlhiv_lt200 = p_hc.prop_lt200(t);
      i_hc.prop_wlhiv_200to350 = 1.0 - p_hc.prop_gte350(t) - p_hc.prop_lt200(t);
      i_hc.prop_wlhiv_gte350 = p_hc.prop_gte350(t);
      i_hc.prop_wlhiv_lt350 = 1 - p_hc.prop_gte350(t);
      i_hc.num_wlhiv = p_hc.mat_hiv_births(t);
    } else {
      i_hc.num_wlhiv_lt200 = 0.0;
      i_hc.num_wlhiv_200to350 = 0.0;
      i_hc.num_wlhiv_gte350 = 0.0;
      i_hc.num_wlhiv = 0.0;
      i_hc.prop_wlhiv_lt200 = 0.0;
      i_hc.prop_wlhiv_200to350 = 0.0;
      i_hc.prop_wlhiv_gte350 = 0.0;
      i_hc.prop_wlhiv_lt350 = 0.0;

      for (int a = 0; a < 35; ++a) {
        i_hc.num_wlhiv_lt200 += n_ha.h_hiv_adult(4, a, 1) + n_ha.h_hiv_adult(5, a, 1) + n_ha.h_hiv_adult(6, a, 1);
        i_hc.num_wlhiv_200to350 += n_ha.h_hiv_adult(3, a, 1) + n_ha.h_hiv_adult(2, a, 1);
        i_hc.num_wlhiv_gte350 += n_ha.h_hiv_adult(0, a, 1) + n_ha.h_hiv_adult(1, a, 1);
      }
      i_hc.num_wlhiv = i_hc.num_wlhiv_200to350 + i_hc.num_wlhiv_gte350 + i_hc.num_wlhiv_lt200;

      if (i_hc.num_wlhiv > 0) {
        i_hc.prop_wlhiv_lt200 = i_hc.num_wlhiv_lt200 / i_hc.num_wlhiv;
        i_hc.prop_wlhiv_200to350 = i_hc.num_wlhiv_200to350 / i_hc.num_wlhiv;
        i_hc.prop_wlhiv_gte350 = i_hc.num_wlhiv_gte350 / i_hc.num_wlhiv;
      } else {
        i_hc.prop_wlhiv_lt200 = 0;
        i_hc.prop_wlhiv_200to350 = 1;
        i_hc.prop_wlhiv_gte350 = 0;
      }
      i_hc.prop_wlhiv_lt350 = i_hc.prop_wlhiv_lt200 + i_hc.prop_wlhiv_200to350;
    }
  };

  void adjust_option_A_B_tr() {
    const auto& p_hc = pars.hc;
    auto& i_hc = intermediate.hc;

    // Option A and B were only authorized for women with greater than 350 CD4, so if the percentage of women
    // on option A/B > the proportion of women in this cd4 category, we assume that some must have a cd4 less than 350
    // option AB will be less effective for these women so we adjust for that
    calc_wlhiv_cd4_proportion();

    auto option_A_B_coverage = i_hc.PMTCT_coverage(0) + i_hc.PMTCT_coverage(1);
    if (option_A_B_coverage > i_hc.prop_wlhiv_gte350) {
      if (i_hc.prop_wlhiv_gte350 > 0) {
        i_hc.excessratio = option_A_B_coverage / i_hc.prop_wlhiv_gte350 - 1;
      } else {
        i_hc.excessratio = 0;
      }
      i_hc.optA_transmission_rate = p_hc.PMTCT_transmission_rate(0, 0, 0) * (1 + i_hc.excessratio);
      i_hc.optB_transmission_rate = p_hc.PMTCT_transmission_rate(0, 1, 0) * (1 + i_hc.excessratio);
    } else {
      i_hc.excessratio = 0.0;
      i_hc.optA_transmission_rate = p_hc.PMTCT_transmission_rate(0, 0, 0) * (1 + i_hc.excessratio);
      i_hc.optB_transmission_rate = p_hc.PMTCT_transmission_rate(0, 1, 0) * (1 + i_hc.excessratio);
    }
  };

  void adjust_option_A_B_bf_tr() {
    const auto& p_hc = pars.hc;
    auto& i_hc = intermediate.hc;

    calc_wlhiv_cd4_proportion();

    // Option A and B were only authorized for women with greater than 350 CD4, so if the percentage of women
    // on option A/B > the proportion of women in this cd4 category, we assume that some must have a cd4 less than 350
    // option AB will be less effective for these women so we adjust for that

    if (i_hc.prop_wlhiv_gte350 > 0) {
      auto option_A_B_coverage = i_hc.PMTCT_coverage(0) + i_hc.PMTCT_coverage(1);
      if (option_A_B_coverage > i_hc.prop_wlhiv_gte350) {
        i_hc.excessratio_bf = option_A_B_coverage - i_hc.prop_wlhiv_gte350;
        auto excess_factor_bf = i_hc.excessratio_bf / option_A_B_coverage * (1.45 / 0.46) +
                                i_hc.prop_wlhiv_gte350;
        i_hc.optA_bf_transmission_rate = excess_factor_bf * p_hc.PMTCT_transmission_rate(4, 0, 1);
        i_hc.optB_bf_transmission_rate = excess_factor_bf * p_hc.PMTCT_transmission_rate(4, 1, 1);
      } else {
        i_hc.optA_bf_transmission_rate = p_hc.PMTCT_transmission_rate(4, 0, 1);
        i_hc.optB_bf_transmission_rate = p_hc.PMTCT_transmission_rate(4, 1, 1);
      }
    } else {
      i_hc.optA_bf_transmission_rate = p_hc.PMTCT_transmission_rate(4, 0, 1);
      i_hc.optB_bf_transmission_rate = p_hc.PMTCT_transmission_rate(4, 1, 1);
    }
  };

  void maternal_incidence_in_pregnancy_tr() {
    const auto& p_dp = pars.dp;
    const auto& p_hc = pars.hc;
    auto& n_dp = state_next.dp;
    auto& n_ha = state_next.ha;
    auto& i_hc = intermediate.hc;

    // Transmission due to incident infections
    i_hc.asfr_sum = 0.0;
    for (int a = 0; a < p_fertility_age_groups; ++a) {
      i_hc.asfr_sum += p_dp.age_specific_fertility_rate(a, t);
    } // end a

    if (p_hc.mat_prev_input(t)) {
      for (int a = 0; a < p_fertility_age_groups; ++a) {
        auto asfr_weight = p_dp.age_specific_fertility_rate(a, t) / i_hc.asfr_sum;
        i_hc.age_weighted_hivneg += asfr_weight * p_hc.adult_female_hivnpop(a, t); // HIV negative 15-49 women weighted for ASFR
        i_hc.age_weighted_infections += asfr_weight * p_hc.adult_female_infections(a, t); // newly infected 15-49 women, weighted for ASFR
      } // end a

      if (i_hc.age_weighted_hivneg > 0.0) {
        i_hc.incidence_rate_wlhiv = i_hc.age_weighted_infections / i_hc.age_weighted_hivneg;
        // 0.75 is 9/12, gestational period, index 7 in the vertical transmission object is the index for maternal seroconversion
        i_hc.perinatal_transmission_from_incidence = i_hc.incidence_rate_wlhiv * 0.75 *
                                                     (p_hc.total_births(t) - i_hc.need_PMTCT) *
                                                     p_hc.vertical_transmission_rate(7, 0);
      } else {
        i_hc.incidence_rate_wlhiv = 0.0;
        i_hc.perinatal_transmission_from_incidence = 0.0;
      }
    } else {
      for (int a = 0; a < p_fertility_age_groups; ++a) {
        auto asfr_weight = p_dp.age_specific_fertility_rate(a, t) / i_hc.asfr_sum;
        i_hc.age_weighted_hivneg += asfr_weight * i_hc.p_hiv_neg_pop(a + 15, 1); // HIV negative 15-49 women weighted for ASFR
        i_hc.age_weighted_infections += asfr_weight * n_ha.p_infections(a + 15, 1); // newly infected 15-49 women, weighted for ASFR
      } // end a

      if (i_hc.age_weighted_hivneg > 0.0) {
        i_hc.incidence_rate_wlhiv = i_hc.age_weighted_infections / i_hc.age_weighted_hivneg;
        //0.75 is 9/12, gestational period, index 7 in the vertical trasnmission object is the index for maternal seroconversion
        i_hc.perinatal_transmission_from_incidence = i_hc.incidence_rate_wlhiv * 0.75 *
                                                     (n_dp.births - i_hc.need_PMTCT) *
                                                     p_hc.vertical_transmission_rate(7, 0);
      } else {
        i_hc.incidence_rate_wlhiv = 0.0;
        i_hc.perinatal_transmission_from_incidence = 0.0;
      }
    }
  };

  void perinatal_tr() {
    const auto& p_hc = pars.hc;
    auto& n_dp = state_next.dp;
    auto& i_hc = intermediate.hc;

    i_hc.births_sum = n_dp.births;

    // TODO: add in patients reallocated
    convert_PMTCT_num_to_perc();
    adjust_option_A_B_tr();
    calc_hiv_negative_pop();

    // Calculate transmission rate
    i_hc.retained_on_ART = i_hc.PMTCT_coverage(4);
    i_hc.retained_started_ART = i_hc.PMTCT_coverage(5);

    // Transmission among women on treatment
    i_hc.perinatal_transmission_rate = i_hc.PMTCT_coverage(0) * i_hc.optA_transmission_rate +
                                       i_hc.PMTCT_coverage(1) * i_hc.optB_transmission_rate +
                                       i_hc.PMTCT_coverage(2) * p_hc.PMTCT_transmission_rate(0, 2, 0) + // SDNVP
                                       i_hc.PMTCT_coverage(3) * p_hc.PMTCT_transmission_rate(0, 3, 0) + //dual ARV
                                       i_hc.retained_on_ART * p_hc.PMTCT_transmission_rate(0, 4, 0) +
                                       i_hc.retained_started_ART * p_hc.PMTCT_transmission_rate(0, 5, 0) +
                                       i_hc.PMTCT_coverage(6) * p_hc.PMTCT_transmission_rate(0, 6, 0);

    i_hc.receiving_PMTCT = i_hc.PMTCT_coverage(0) + i_hc.PMTCT_coverage(1) +
                           i_hc.PMTCT_coverage(2) + i_hc.PMTCT_coverage(3) +
                           i_hc.retained_on_ART + i_hc.retained_started_ART +
                           i_hc.PMTCT_coverage(6);

    i_hc.no_PMTCT = 1 - i_hc.receiving_PMTCT;
    i_hc.no_PMTCT = std::max(i_hc.no_PMTCT, 0.0);

    // Transmission among women not on treatment
    if (i_hc.num_wlhiv > 0) {
      auto untreated_vertical_tr = i_hc.prop_wlhiv_lt200 * p_hc.vertical_transmission_rate(4, 0) +
                                   i_hc.prop_wlhiv_200to350 * p_hc.vertical_transmission_rate(2, 0) +
                                   i_hc.prop_wlhiv_gte350 * p_hc.vertical_transmission_rate(0, 0);
      i_hc.perinatal_transmission_rate += i_hc.no_PMTCT * untreated_vertical_tr;
    }
    i_hc.perinatal_transmission_rate_bf_calc = i_hc.perinatal_transmission_rate;

    maternal_incidence_in_pregnancy_tr();

    if (i_hc.need_PMTCT > 0.0) {
      i_hc.perinatal_transmission_rate += i_hc.perinatal_transmission_from_incidence / i_hc.need_PMTCT;
    }
  };

  void maternal_incidence_in_bf_tr() {
    const auto& p_hc = pars.hc;
    auto& i_hc = intermediate.hc;

    for (int bf = 0; bf < hBF; ++bf) {
      i_hc.bf_at_risk += i_hc.incidence_rate_wlhiv / 12 * 2 *
                         (1 - p_hc.breastfeeding_duration_no_art(bf, t));
    }
    i_hc.bf_incident_hiv_transmission_rate = i_hc.bf_at_risk * p_hc.vertical_transmission_rate(7, 1);
  };

  void bf_dropout(int bf) {
    const auto& p_hc = pars.hc;
    auto& i_hc = intermediate.hc;

    const auto hPS_dropout_idx = (bf < 6) ? 4 : 5;
    const auto PMTCT_retention = 1 - p_hc.PMTCT_dropout(hPS_dropout_idx, t) * 2;

    i_hc.PMTCT_coverage(0) *= PMTCT_retention; //opt A
    i_hc.PMTCT_coverage(1) *= PMTCT_retention; //opt B
    i_hc.PMTCT_coverage(4) *= PMTCT_retention; //before pregnancy
    i_hc.PMTCT_coverage(5) *= PMTCT_retention; //>4 weeks
    i_hc.PMTCT_coverage(6) *= PMTCT_retention; //<4 weeks
  };

  void run_bf_transmission_rate(int bf_start, int bf_end, int index) {
    const auto& p_hc = pars.hc;
    auto& i_hc = intermediate.hc;

    for (int bf = bf_start; bf < bf_end; bf++) {
      if (bf == 0) {
        // Perinatal transmission accounts for transmission up to 6 weeks, so we only use 1/4 of
        // transmission from the first breastfeeding period
        i_hc.bf_scalar = 0.25;
      } else {
        i_hc.bf_scalar = 1.0;
        // dropout only occurs after the first month of breastfeeding
        bf_dropout(bf);
      }

      // i_hc.perinatal_transmission_rate_bf_calc is the transmission that has already occurred due to perinatal transmission
      // i_hc.percent_no_treatment is the percentage of women who are still vulnerable to HIV transmission to their babies
      i_hc.percent_no_treatment = 1 - i_hc.perinatal_transmission_rate_bf_calc - i_hc.bf_transmission_rate(index);

      if (index > 0) {
        for (int bf = 0; bf < index; ++bf) {
          i_hc.percent_no_treatment -= i_hc.bf_transmission_rate(bf);
        }
      }

      for (int hp = 0; hp < hPS; hp++) {
        i_hc.percent_on_treatment = 0;
        i_hc.percent_no_treatment -=  i_hc.PMTCT_coverage(hp);

        if (hp <= 1) continue;

        // sdnvp stratifies transmission by CD4, but spectrum only uses one
        const auto hDS_idx = (hp == 2) ? 0 : 4;
        auto tr = i_hc.PMTCT_coverage(hp) * p_hc.PMTCT_transmission_rate(hDS_idx, hp, 1) *
                  2 * (1 - p_hc.breastfeeding_duration_art(bf, t)) * i_hc.bf_scalar;
        i_hc.PMTCT_coverage(hp) -= tr;
        i_hc.bf_transmission_rate(index) += tr;
      }

      // No treatment
      if (p_hc.breastfeeding_duration_no_art(bf, t) < 1) {
        i_hc.percent_no_treatment = std::max(i_hc.percent_no_treatment, 0.0);
        auto untreated_vertical_bf_tr = i_hc.prop_wlhiv_lt200 * p_hc.vertical_transmission_rate(4, 1) +
                                        i_hc.prop_wlhiv_200to350 * p_hc.vertical_transmission_rate(2, 1) +
                                        i_hc.prop_wlhiv_gte350 * p_hc.vertical_transmission_rate(0, 1);
        i_hc.bf_transmission_rate(index) += i_hc.bf_scalar * i_hc.percent_no_treatment *
                                            untreated_vertical_bf_tr *
                                            2 * (1 - p_hc.breastfeeding_duration_no_art(bf, t));
      }
    }
  };

  void nosocomial_infections() {
    const auto& p_hc = pars.hc;
    auto& n_ha = state_next.ha;
    auto& n_hc = state_next.hc;

    for (int s = 0; s < NS; ++s) {
      // Run only first 5 age groups in total population 0, 1, 2, 3, 4
      for (int a = 0; a < hc2_agestart; ++a) {
        if (p_hc.hc_nosocomial(t) > 0) {
          // 5.0 is used because we want to evenly distribute across the 5 age groups in 0-4
          n_ha.p_infections(a, s) = p_hc.hc_nosocomial(t) / (5.0 * NS);
          // Putting all nosocomial acquired HIV infections in perinatally acquired infection timing and highest CD4 category to match Spectrum implementation
          n_hc.hc1_hiv_pop(0, 0, a, s) += n_ha.p_infections(a, s);
        }
      } // end a
    } // end NS
  };

  void add_infections() {
    const auto& p_dp = pars.dp;
    const auto& p_hc = pars.hc;
    auto& n_ha = state_next.ha;
    auto& n_dp = state_next.dp;
    auto& n_hc = state_next.hc;
    auto& i_hc = intermediate.hc;

    if (n_hc.hiv_births > 0) {
      perinatal_tr();

      // Perinatal transmission
      auto perinatal_transmission_births = n_hc.hiv_births * i_hc.perinatal_transmission_rate;
      for (int s = 0; s < NS; ++s) {
        for (int hd = 0; hd < hc1DS; ++hd) {
          n_hc.hc1_hiv_pop(hd, 0, 0, s) += perinatal_transmission_births *
                                           p_dp.births_sex_prop(s, t) * p_hc.hc1_cd4_dist(hd);
        } // end hc1DS
        auto perinatal_births_by_sex = perinatal_transmission_births * p_dp.births_sex_prop(s, t);
        n_ha.p_infections(0, s) += perinatal_births_by_sex;
        n_hc.infection_by_type(0, 0, s) += perinatal_births_by_sex;
      } // end NS

      // Breastfeeding transmission
      // 0-6
      maternal_incidence_in_bf_tr();
      adjust_option_A_B_bf_tr();
      convert_PMTCT_pre_bf();
      run_bf_transmission_rate(0, 3, 0);

      real_type total_births = 0.0;
      if (p_hc.mat_prev_input(t)) {
        total_births = p_hc.total_births(t);
      } else {
        total_births = n_dp.births;
      }

      // 0-6
      for (int s = 0; s < NS; ++s) {
        auto bf_hiv_by_sex = n_hc.hiv_births * p_dp.births_sex_prop(s, t) *
                             i_hc.bf_transmission_rate(0);
        // vertical infection from maternal infection during breastfeeding
        bf_hiv_by_sex += (total_births - n_hc.hiv_births) * p_dp.births_sex_prop(s, t) *
                         i_hc.bf_incident_hiv_transmission_rate;
        for (int hd = 0; hd < hc1DS; ++hd) {
          n_hc.hc1_hiv_pop(hd, 1, 0, s) += p_hc.hc1_cd4_dist(hd) * bf_hiv_by_sex;
        } // end hc1DS
        n_ha.p_infections(0, s) += bf_hiv_by_sex;
        n_hc.infection_by_type(1, 0, s) += bf_hiv_by_sex;
      } // end NS

      // 6-12
      run_bf_transmission_rate(3, 6, 1);
      for (int s = 0; s < NS; ++s) {
        auto bf_hiv_by_sex = n_hc.hiv_births * p_dp.births_sex_prop(s, t) * i_hc.bf_transmission_rate(1);
        for (int hd = 0; hd < hc1DS; ++hd) {
          n_hc.hc1_hiv_pop(hd, 2, 0, s) += p_hc.hc1_cd4_dist(hd) * bf_hiv_by_sex;
        } // end hc1DS
        n_ha.p_infections(0, s) += bf_hiv_by_sex;
        n_hc.infection_by_type(2, 0, s) += bf_hiv_by_sex;
      } // end NS

      // 12 plus
      run_bf_transmission_rate(6, 12, 2);
      run_bf_transmission_rate(12, hBF, 3);
      for (int s = 0; s < NS; ++s) {
        auto total_pop_12_24 = n_dp.p_total_pop(1, 0) - n_ha.p_hiv_pop(1, 0) +
                               n_dp.p_total_pop(1, 1) - n_ha.p_hiv_pop(1, 1);
        if (s == 1) {
          total_pop_12_24 -= n_ha.p_infections(1, 0);
        }
        auto uninfected_prop_12_24 = (n_dp.p_total_pop(1, s) - n_ha.p_hiv_pop(1, s)) / total_pop_12_24;

        auto total_pop_24_plus = n_dp.p_total_pop(2, 0) - n_ha.p_hiv_pop(2, 0) +
                                 n_dp.p_total_pop(2, 1) - n_ha.p_hiv_pop(2, 1);
        if (s == 1) {
          total_pop_24_plus -= n_ha.p_infections(2, 0);
        }
        auto uninfected_prop_24_plus = (n_dp.p_total_pop(2, s) - n_ha.p_hiv_pop(2, s)) / total_pop_24_plus;

        for (int hd = 0; hd < hc1DS; ++hd) {
          auto bf_hiv_transmission_12_24 = n_hc.hiv_births * i_hc.bf_transmission_rate(2) * uninfected_prop_12_24;
          auto bf_hiv_transmission_24_plus = n_hc.hiv_births * i_hc.bf_transmission_rate(3) * uninfected_prop_24_plus;
          // 12-24
          n_hc.hc1_hiv_pop(hd, 3, 1, s) += p_hc.hc1_cd4_dist(hd) * bf_hiv_transmission_12_24;
          n_ha.p_infections(1, s) += p_hc.hc1_cd4_dist(hd) * bf_hiv_transmission_12_24;
          n_hc.infection_by_type(3, 1, s) += p_hc.hc1_cd4_dist(hd) * bf_hiv_transmission_12_24;

          // 24 plus
          n_hc.hc1_hiv_pop(hd, 3, 2, s) += p_hc.hc1_cd4_dist(hd) * bf_hiv_transmission_24_plus;
          n_ha.p_infections(2, s) += p_hc.hc1_cd4_dist(hd) * bf_hiv_transmission_24_plus;
          n_hc.infection_by_type(3, 2, s) += p_hc.hc1_cd4_dist(hd) * bf_hiv_transmission_24_plus;
        } // end hc1DS
      } // end NS
    }
  };

  void art_eligibility_by_age() {
    const auto& p_hc = pars.hc;
    auto& n_hc = state_next.hc;

    // all children under a certain age eligible for ART
    for (int s = 0; s < NS; ++s) {
      for (int a = 0; a < p_hc.hc_art_elig_age(t); ++a) {
        for (int cat = 0; cat < hcTT; ++cat) {
          for (int hd = 0; hd < hc1DS; ++hd) {
            if (a < hc2_agestart) {
              n_hc.hc_art_need_init(hd, cat, a, s) += n_hc.hc1_hiv_pop(hd, cat, a, s);
            } else if (hd < hc2DS) {
              n_hc.hc_art_need_init(hd, cat, a, s) += n_hc.hc2_hiv_pop(hd, cat, a - hc2_agestart, s);
            }
          } // end hc1DS
        } // end a
      } // end hcTT
    } // end NS
  };

  void art_eligibility_by_cd4() {
    const auto& p_hc = pars.hc;
    auto& n_hc = state_next.hc;

    // all children under a certain CD4 eligible for ART
    for (int s = 0; s < NS; ++s) {
      for (int cat = 0; cat < hcTT; ++cat) {
        for (int a = p_hc.hc_art_elig_age(t); a < p_idx_hiv_first_adult; ++a) {
          for (int hd = 0; hd < hc1DS; ++hd) {
            if (hd >= p_hc.hc_art_elig_cd4(a, t)) {
              if (a < hc2_agestart) {
                n_hc.hc_art_need_init(hd, cat, a, s) += n_hc.hc1_hiv_pop(hd, cat, a, s);
              } else if (hd < hc2DS) {
                n_hc.hc_art_need_init(hd, cat, a, s) += n_hc.hc2_hiv_pop(hd, cat, a - hc2_agestart, s);
              }
            }
          } // end hc1DS
        } // end a
      } // end hcTT
    } // end NS
  };

  void need_for_cotrim() {
    const auto& p_hc = pars.hc;
    const auto& c_hc = state_curr.hc;
    auto& n_hc = state_next.hc;

    // Births from the last 18 months are eligible
    n_hc.ctx_need = n_hc.hiv_births * 1.5;

    // All children 1.5-4 eligible
    for (int s = 0; s < NS; ++s) {
      for (int a = 1; a < hc2_agestart; ++a) {
        for (int cat = 0; cat < hcTT; ++cat) {
          for (int hd = 0; hd < hc1DS; ++hd) {
            if (a == 1) {
              n_hc.ctx_need += n_hc.hc1_hiv_pop(hd, cat, a, s) * 0.5;
            } else {
              n_hc.ctx_need += n_hc.hc1_hiv_pop(hd, cat, a, s);
            }
          }
        }
      }
    }

    // Children under five on ART also eligible
    for (int s = 0; s < NS; ++s) {
      for (int a = 0; a < hc2_agestart; ++a) {
        for (int dur = 0; dur < hTS; ++dur) {
          for (int hd = 0; hd < hc1DS; ++hd) {
            if (a == 1) {
              n_hc.ctx_need += n_hc.hc1_art_pop(dur, hd, a, s) * 0.5;
            } else {
              n_hc.ctx_need += n_hc.hc1_art_pop(dur, hd, a, s) ;
            }
          } // end hc1DS
        } // end a
      } // end hcTT
    } // end NS

    // All ART eligible children ages 5-14 eligible
    // Spectrum uses a lagged population and eligibility for children over five (TODO: verify)
    for (int s = 0; s < NS; ++s) {
      for (int cat = 0; cat < hcTT; ++cat) {
        for (int a = hc2_agestart; a < p_idx_hiv_first_adult; ++a) {
          for (int hd = 0; hd < hc2DS; ++hd) {
            if (a < p_hc.hc_art_elig_age(t) || hd >= p_hc.hc_art_elig_cd4(a, t - 1)) {
              n_hc.ctx_need += c_hc.hc2_hiv_pop(hd, cat, a - hc2_agestart, s);
            }
          } // end hc1DS
        } // end a
      } // end hcTT
    } // end NS
  };

  void get_cotrim_effect(int art_flag) {
    const auto& p_hc = pars.hc;
    auto& n_hc = state_next.hc;
    auto& i_hc = intermediate.hc;

    if (p_hc.ctx_val_is_percent(t)) {
      i_hc.ctx_mean(art_flag) = 1 - p_hc.ctx_effect(art_flag) * p_hc.ctx_val(t);
    } else if (n_hc.ctx_need > 0) {
      i_hc.ctx_mean(art_flag) = 1 - p_hc.ctx_effect(art_flag) * p_hc.ctx_val(t) / n_hc.ctx_need;
    }
  };

  void cd4_mortality() {
    const auto& p_hc = pars.hc;
    auto& n_hc = state_next.hc;
    auto& i_hc = intermediate.hc;
    auto art_flag = 0;

    get_cotrim_effect(art_flag);

    for (int s = 0; s < NS; ++s) {
      for (int a = 0; a < hc2_agestart; ++a) {
        for (int cat = 0; cat < hcTT; ++cat) {
          for (int hd = 0; hd < hc1DS; ++hd) {
            auto hiv_deaths_strat = i_hc.ctx_mean(art_flag) * n_hc.hc1_hiv_pop(hd, cat, a, s) * p_hc.hc1_cd4_mort(hd, cat, a);
            i_hc.hc_posthivmort(hd, cat, a, s) = n_hc.hc1_hiv_pop(hd, cat, a, s) - hiv_deaths_strat;
          }
        }
      }
    }

    for (int s = 0; s < NS; ++s) {
      for (int a = hc2_agestart; a < p_idx_hiv_first_adult; ++a) {
        for (int cat = 0; cat < hcTT; ++cat) {
          for (int hd = 0; hd < hc2DS; ++hd) {
            auto hiv_deaths_strat = i_hc.ctx_mean(art_flag) * n_hc.hc2_hiv_pop(hd, cat, a - hc2_agestart, s) *
                                    p_hc.hc2_cd4_mort(hd, cat, a - hc2_agestart);
            i_hc.hc_posthivmort(hd, cat, a, s) = n_hc.hc2_hiv_pop(hd, cat, a - hc2_agestart, s) -
                                                hiv_deaths_strat;
          }
        }
      }
    }

    // progress through CD4 categories
    for (int s = 0; s < NS; ++s) {
      for (int hd = 1; hd < hc1DS; ++hd) {
        for (int a = 0; a < hc2_agestart; ++a) {
          for (int cat = 0; cat < hcTT; ++cat) {
            const auto& coarse_hc1_cd4_prog = p_hc.hc1_cd4_prog(hd - 1, hc_age_coarse_cd4[a], s);
            auto cd4_grad = coarse_hc1_cd4_prog *
                            (i_hc.hc_posthivmort(hd - 1, cat, a, s) + n_hc.hc1_hiv_pop(hd - 1, cat, a, s)) /
                            2.0;
            i_hc.hc_grad(hd - 1, cat, a, s) -= cd4_grad; // moving to next cd4 category
            i_hc.hc_grad(hd, cat, a, s) += cd4_grad; // moving into this cd4 category
          }
        }
      }
    }

    // progress through CD4 categories
    for (int s = 0; s < NS; ++s) {
      for (int hd = 1; hd < hc2DS; ++hd) {
        for (int a = hc2_agestart; a < p_idx_hiv_first_adult; ++a) {
          for (int cat = 0; cat < hcTT; ++cat) {
            auto cd4_grad = p_hc.hc2_cd4_prog(hd - 1, 0, s) *
                            (i_hc.hc_posthivmort(hd - 1, cat, a, s) + n_hc.hc2_hiv_pop(hd - 1, cat, a - hc2_agestart, s)) /
                            2.0;
            i_hc.hc_grad(hd - 1, cat, a, s) -= cd4_grad; // moving to next cd4 category
            i_hc.hc_grad(hd, cat, a, s) += cd4_grad; // moving into this cd4 category
          }
        }
      }
    }
  };

  void run_child_hiv_mort() {
    const auto& p_hc = pars.hc;
    auto& n_hc = state_next.hc;
    auto& i_hc = intermediate.hc;
    auto art_flag = 0;

    for (int s = 0; s < NS; ++s) {
      for (int a = 0; a < hc2_agestart; ++a) {
        for (int cat = 0; cat < hcTT; ++cat) {
          for (int hd = 0; hd < hc1DS; ++hd) {
            auto cd4_deaths_grad = i_hc.ctx_mean(art_flag) * n_hc.hc1_hiv_pop(hd, cat, a, s) *
                                   p_hc.hc1_cd4_mort(hd, cat, a);
            i_hc.hc_grad(hd, cat, a, s) -= cd4_deaths_grad;
            n_hc.hc1_noart_aids_deaths(hd, cat, a, s) += cd4_deaths_grad;
          }
        }
      }
    }

    for (int s = 0; s < NS; ++s) {
      for (int a = hc2_agestart; a < p_idx_hiv_first_adult; ++a) {
        for (int cat = 0; cat < hcTT; ++cat) {
          for (int hd = 0; hd < hc2DS; ++hd) {
            auto cd4_mort_grad = i_hc.ctx_mean(art_flag) *
                                 n_hc.hc2_hiv_pop(hd, cat, a - hc2_agestart, s) *
                                 p_hc.hc2_cd4_mort(hd, cat, a - hc2_agestart);
            i_hc.hc_grad(hd, cat, a, s) -= cd4_mort_grad;
            n_hc.hc2_noart_aids_deaths(hd, cat, a - hc2_agestart, s) += cd4_mort_grad;
          }
        }
      }
    }
  };

  void add_child_grad() {
    auto& n_hc = state_next.hc;
    auto& i_hc = intermediate.hc;

    // add on transitions
    for (int s = 0; s < NS; ++s) {
      for (int a = 0; a < hc2_agestart; ++a) {
        for (int cat = 0; cat < hcTT; ++cat) {
          for (int hd = 0; hd < hc1DS; ++hd) {
            n_hc.hc1_hiv_pop(hd, cat, a, s) += i_hc.hc_grad(hd, cat, a, s);
          } // end hc1DS
        } // end cat
      } // end a
    } // end s

    for (int s = 0; s < NS; ++s) {
      for (int a = hc2_agestart; a < p_idx_hiv_first_adult; ++a) {
        for (int cat = 0; cat < hcTT; ++cat) {
          for (int hd = 0; hd < hc2DS; ++hd) {
            n_hc.hc2_hiv_pop(hd, cat, a - hc2_agestart, s) += i_hc.hc_grad(hd, cat, a, s);
          } //end hc2DS
        } //end cat
      } //end a
    } // end s
  };

  void eligible_for_treatment() {
    auto& n_hc = state_next.hc;
    auto& i_hc = intermediate.hc;

    art_eligibility_by_age();
    art_eligibility_by_cd4();

    for (int s = 0; s < NS; ++s) {
      for (int cat = 0; cat < hcTT; ++cat) {
        for (int a = 0; a < p_idx_hiv_first_adult; ++a) {
          for (int hd = 0; hd < hc1DS; ++hd) {
            i_hc.eligible(hd, a, s) += n_hc.hc_art_need_init(hd, cat, a, s);
          } // end hc1DS
        } // end a
      } // end hcTT
    } // end NS
  };

  void on_art_mortality(int t_art_idx, int art_flag) {
    const auto& p_hc = pars.hc;
    auto& n_hc = state_next.hc;
    auto& i_hc = intermediate.hc;

    get_cotrim_effect(art_flag);

    for (int s = 0; s < NS; ++s) {
      for (int a = 0; a < p_idx_hiv_first_adult; ++a) {
        for (int hd = 0; hd < hc1DS; ++hd) {
          i_hc.hc_death_rate = 0.0;
          i_hc.hc_art_grad(t_art_idx, hd, a, s) = 0.0;
          if (t_art_idx == 0) {
            if (a < hc2_agestart) {
              i_hc.hc_death_rate = p_hc.hc_art_mort_rr(t_art_idx, a, t) *
                                   (p_hc.hc1_art_mort(hd, 0, a) + p_hc.hc1_art_mort(hd, 1, a)) /
                                   2.0;
            } else if (hd < hc2DS) {
              bool is_hc2_art_pop = n_hc.hc2_art_pop(t_art_idx, hd, a - hc2_agestart, s) > 0;
              if (is_hc2_art_pop) {
                i_hc.hc_death_rate = p_hc.hc_art_mort_rr(t_art_idx, a, t) *
                                     (p_hc.hc2_art_mort(hd, 0, a - hc2_agestart) + p_hc.hc2_art_mort(hd, 1, a - hc2_agestart)) /
                                     2.0;
              }
            }
          } else {
            if (a < hc2_agestart) {
              i_hc.hc_death_rate = p_hc.hc_art_mort_rr(t_art_idx, a, t) * p_hc.hc1_art_mort(hd, 2, a);
            } else if (hd < hc2DS) {
              bool is_hc2_art_pop = n_hc.hc2_art_pop(t_art_idx, hd, a - hc2_agestart, s) > 0;
              if (is_hc2_art_pop) {
                i_hc.hc_death_rate = p_hc.hc_art_mort_rr(t_art_idx, a, t) * p_hc.hc2_art_mort(hd, 2, a - hc2_agestart);
              }
            }
          }

          // ctx reduction on mortality for those on ART
          i_hc.hc_death_rate *= i_hc.ctx_mean(art_flag);

          if (a < hc2_agestart) {
            bool any_hc1_art_deaths = i_hc.hc_death_rate * n_hc.hc1_art_pop(t_art_idx, hd, a, s) >= 0;
            if (any_hc1_art_deaths) {
              i_hc.hc_art_grad(t_art_idx, hd, a, s) -= i_hc.hc_death_rate * n_hc.hc1_art_pop(t_art_idx, hd, a, s);
              n_hc.hc1_art_pop(t_art_idx, hd, a, s) += i_hc.hc_art_grad(t_art_idx, hd, a, s);
              n_hc.hc1_art_aids_deaths(t_art_idx, hd, a, s) -= i_hc.hc_art_grad(t_art_idx, hd, a, s);
            }
          } else if (hd < hc2DS) {
            bool any_hc2_art_deaths = i_hc.hc_death_rate * n_hc.hc2_art_pop(t_art_idx, hd, a - hc2_agestart, s) >= 0;
            if (any_hc2_art_deaths) {
              i_hc.hc_art_grad(t_art_idx, hd, a, s) -= i_hc.hc_death_rate *
                                                       n_hc.hc2_art_pop(t_art_idx, hd, a - hc2_agestart, s);
              n_hc.hc2_art_pop(t_art_idx, hd, a - hc2_agestart, s) += i_hc.hc_art_grad(t_art_idx, hd, a, s);
              n_hc.hc2_art_aids_deaths(t_art_idx, hd, a - hc2_agestart, s) -= i_hc.hc_art_grad(t_art_idx, hd, a, s);
            }
          }
        } // end a
      } // end hc1DS
    } // end NS
  };

  void deaths_this_year() {
    const auto& p_hc = pars.hc;
    auto& n_hc = state_next.hc;
    auto& i_hc = intermediate.hc;

    for (int dur = 0; dur < hTS; ++dur) {
      for (int s = 0; s < NS; ++s) {
        for (int hd = 0; hd < hc1DS; ++hd) {
          for (int a = 0; a < p_idx_hiv_first_adult; ++a) {
            i_hc.hc_death_rate = 0.0;

            if (dur == 0) {
              if (a < hc2_agestart) {
                i_hc.hc_death_rate = p_hc.hc_art_mort_rr(dur, a, t) *
                                     (p_hc.hc1_art_mort(hd, 0, a) + p_hc.hc1_art_mort(hd, 1, a)) /
                                     2.0;
              } else if (hd < hc2DS) {
                bool is_hc2_art_pop = n_hc.hc2_art_pop(dur, hd, a - hc2_agestart, s) > 0;
                if (is_hc2_art_pop) {
                  i_hc.hc_death_rate = p_hc.hc_art_mort_rr(dur, a, t) *
                                       (p_hc.hc2_art_mort(hd, 0, a - hc2_agestart) + p_hc.hc2_art_mort(hd, 1, a - hc2_agestart)) /
                                       2.0;
                }
              }
            } else {
              if (a < hc2_agestart) {
                i_hc.hc_death_rate = p_hc.hc_art_mort_rr(dur, a, t) * p_hc.hc1_art_mort(hd, 2, a);
              } else if (hd < hc2DS) {
                bool is_hc2_art_pop = n_hc.hc2_art_pop(dur, hd, a - hc2_agestart, s) > 0;
                if (is_hc2_art_pop) {
                  i_hc.hc_death_rate = p_hc.hc_art_mort_rr(dur, a, t) * p_hc.hc2_art_mort(hd, 2, a - hc2_agestart);
                }
              }
            }

            // ctx reduction on mortality for those on ART
            // NOTE: ART initiation calculations don't include the effect of cotrim (TODO: verify)
            // i_hc.hc_death_rate *= i_hc.ctx_mean(art_flag);
            if (a < hc2_agestart) {
              bool any_hc1_art_deaths = i_hc.hc_death_rate * n_hc.hc1_art_pop(dur, hd, a, s) >= 0;
              if (any_hc1_art_deaths) {
                i_hc.hc_art_deaths(hc_age_coarse[a]) += i_hc.hc_death_rate * n_hc.hc1_art_pop(dur, hd, a, s);
              }
            } else if (hd < hc2DS) {
              bool any_hc2_art_deaths = i_hc.hc_death_rate * n_hc.hc2_art_pop(dur, hd, a - hc2_agestart, s) >= 0;
              if (any_hc2_art_deaths) {
                i_hc.hc_art_deaths(hc_age_coarse[a]) += i_hc.hc_death_rate *
                                                             n_hc.hc2_art_pop(dur, hd, a - hc2_agestart, s);
              }
            }
          } // end a
        } // end hc1DS
      } // end NS
    } // end dur
    i_hc.hc_art_deaths(0) = i_hc.hc_art_deaths(1) + i_hc.hc_art_deaths(2) + i_hc.hc_art_deaths(3);
  };

  void progress_time_on_art(int curr_t_idx, int end_t_idx) {
    auto& n_hc = state_next.hc;

    // Progress ART to the correct time on ART
    for (int hd = 0; hd < hc1DS; ++hd) {
      for (int a = 0; a < p_idx_hiv_first_adult; ++a) {
        for (int s = 0; s < NS; ++s) {
          if (a < hc2_agestart) {
            if (n_hc.hc1_art_pop(curr_t_idx, hd, a, s) > 0) {
              n_hc.hc1_art_pop(end_t_idx, hd, a, s) += n_hc.hc1_art_pop(curr_t_idx, hd, a, s);
              n_hc.hc1_art_pop(curr_t_idx, hd, a, s) -= n_hc.hc1_art_pop(curr_t_idx, hd, a, s);
            }
          } else if (hd < hc2DS) {
            n_hc.hc2_art_pop(end_t_idx, hd, a - hc2_agestart, s) += n_hc.hc2_art_pop(curr_t_idx, hd, a - hc2_agestart, s);
            n_hc.hc2_art_pop(curr_t_idx, hd, a - hc2_agestart, s) -= n_hc.hc2_art_pop(curr_t_idx, hd, a - hc2_agestart, s);
          }
        } // end NS
      } // end a
    } // end hc1DS
  };

  void calc_on_art() {
    auto& n_hc = state_next.hc;
    auto& i_hc = intermediate.hc;

    for (int s = 0; s < NS; ++s) {
      for (int a = 0; a < p_idx_hiv_first_adult; ++a) {
        for (int hd = 0; hd < hc1DS; ++hd) {
          for (int dur = 0; dur < hTS; ++dur) {
            if (a < hc2_agestart) {
              i_hc.on_art(hc_age_coarse[a]) += n_hc.hc1_art_pop(dur, hd, a, s);
            } else if (hd < (hc2DS)) {
              i_hc.on_art(hc_age_coarse[a]) += n_hc.hc2_art_pop(dur, hd, a - hc2_agestart, s);
            }
          }
        }
      }
    }
  };

  void calc_total_and_unmet_art_need() {
    auto& i_hc = intermediate.hc;

    for (int s = 0; s < NS; ++s) {
      for (int a = 0; a < p_idx_hiv_first_adult; ++a) {
        for (int hd = 0; hd < hc1DS; ++hd) {
          i_hc.unmet_need(hc_age_coarse[a]) += i_hc.eligible(hd, a, s);
        } // end hc1DS
      } // end a
    } // end NS

    for (int ag = 1; ag < hcAG_coarse; ++ag) {
      i_hc.on_art(0) += i_hc.on_art(ag);
      i_hc.unmet_need(0) += i_hc.unmet_need(ag);
      i_hc.total_need(ag) += i_hc.on_art(ag) + i_hc.unmet_need(ag) + i_hc.hc_art_deaths(ag);
    } // end ag
    i_hc.total_need(0) = i_hc.on_art(0) + i_hc.unmet_need(0) + i_hc.hc_art_deaths(0);
  };

  void age_specific_art_last_year() {
    const auto& p_hc = pars.hc;
    const auto& c_hc = state_curr.hc;
    auto& i_hc = intermediate.hc;

    if (p_hc.hc_art_is_age_spec(t - 1)) {
      for (int ag = 1; ag < hcAG_coarse; ++ag) {
        i_hc.total_art_last_year(ag) = p_hc.hc_art_val(ag, t - 1);
      } // end ag
    } else {
      for (int s = 0; s < NS; ++s) {
        for (int a = 0; a < p_idx_hiv_first_adult; ++a) {
          for (int hd = 0; hd < hc1DS; ++hd) {
            for (int dur = 0; dur < hTS; ++dur) {
              if (a < hc2_agestart) {
                i_hc.total_art_last_year(hc_age_coarse[a]) += c_hc.hc1_art_pop(dur, hd, a, s);
              } else if (hd < hc2DS) {
                i_hc.total_art_last_year(hc_age_coarse[a]) += c_hc.hc2_art_pop(dur, hd, a - hc2_agestart, s);
              }
            }
          }
        }
      }
      i_hc.total_art_last_year(0) = i_hc.total_art_last_year(1) +
                                    i_hc.total_art_last_year(2) +
                                    i_hc.total_art_last_year(3);

      for (int ag = 1; ag < hcAG_coarse; ++ag) {
        i_hc.total_art_last_year(ag) = p_hc.hc_art_val(0, t - 1) *
                                      i_hc.total_art_last_year(ag) / i_hc.total_art_last_year(0);
        if (p_hc.hc_art_isperc(t - 1)) {
          i_hc.total_art_last_year(ag) *= i_hc.total_need(0) + i_hc.hc_art_deaths(ag);
        }
      } // end ag
    }
  };

  void art_last_year() {
    const auto& p_hc = pars.hc;
    auto& i_hc = intermediate.hc;

    if (p_hc.hc_art_is_age_spec(t)) {
      // If the present time step is age specific, we need to calculate what last years age spec
      // breakdown would have been

      // Age specific ART will always be entered as a number
      age_specific_art_last_year();
    } else {
      if (p_hc.hc_art_isperc(t - 1)) {
        // ART entered as percent last year so convert to number
        i_hc.total_art_last_year(0) = p_hc.hc_art_val(0, t - 1) * i_hc.total_need(0);
      } else if (p_hc.hc_art_is_age_spec(t - 1)) {
        // ART entered as number last year but this year isn't then aggregate
        // ages
        i_hc.total_art_last_year(0) = p_hc.hc_art_val(1, t - 1) +
                                      p_hc.hc_art_val(2, t - 1) +
                                      p_hc.hc_art_val(3, t - 1);
      } else {
        // Last year was age aggregated and a number so use previous value
        i_hc.total_art_last_year(0) = p_hc.hc_art_val(0, t - 1);
      }
    }
  };

  void art_this_year() {
    const auto& p_hc = pars.hc;
    auto& i_hc = intermediate.hc;

    for (int ag = 0; ag < hcAG_coarse; ++ag) {
      i_hc.total_art_this_year(ag) = p_hc.hc_art_val(ag, t);
      if (p_hc.hc_art_isperc(t)) {
        i_hc.total_art_this_year(ag) *= i_hc.total_need(ag);
      }
    } // end ag
  };

  void calc_art_initiates() {
    const auto& p_hc = pars.hc;
    auto& n_hc = state_next.hc;
    auto& i_hc = intermediate.hc;

    calc_on_art();
    deaths_this_year();
    calc_total_and_unmet_art_need();
    art_last_year();
    art_this_year();

    i_hc.retained = 1 - p_hc.hc_art_ltfu(t);
    for (int ag = 0; ag < hcAG_coarse; ++ag) {
      auto average_art_by_year = (i_hc.total_art_last_year(ag) + i_hc.total_art_this_year(ag)) /
                                2.0;
      n_hc.hc_art_init(ag) = std::max(i_hc.hc_art_deaths(ag) + average_art_by_year - i_hc.on_art(ag) * i_hc.retained, 0.0);
      n_hc.hc_art_init(ag) = std::min(n_hc.hc_art_init(ag),
                                      i_hc.unmet_need(ag) + i_hc.on_art(ag) * p_hc.hc_art_ltfu(t));
    } // end ag
  };

  void art_ltfu() {
    const auto& p_hc = pars.hc;
    auto& n_hc = state_next.hc;
    auto& i_hc = intermediate.hc;

    for (int s = 0; s < NS; ++s) {
      for (int a = 0; a < p_idx_hiv_first_adult; ++a) {
        for (int hd = 0; hd < hc1DS; ++hd) {
          for (int cat = 0; cat < hcTT; ++cat) {
            if (a < hc2_agestart) {
              i_hc.hc_hiv_total(hd, a, s) += n_hc.hc1_hiv_pop(hd, cat, a, s);
            } else if (hd < hc2DS) {
              i_hc.hc_hiv_total(hd, a, s) += n_hc.hc2_hiv_pop(hd, cat, a - hc2_agestart, s);
            }
          } // end hcTT
        } // end hc1DS
      } // end a
    } // end NS

    for (int s = 0; s <NS; ++s) {
      for (int a = 0; a < p_idx_hiv_first_adult; ++a) {
        for (int hd = 0; hd < hc1DS; ++hd) {
          for (int cat = 0; cat < hcTT; ++cat) {
            if (a < hc2_agestart) {
              i_hc.hc_hiv_dist(hd, cat, a, s) += n_hc.hc1_hiv_pop(hd, cat, a, s) / i_hc.hc_hiv_total(hd, a, s);
            } else if (hd < hc2DS) {
              i_hc.hc_hiv_dist(hd, cat, a, s) += n_hc.hc2_hiv_pop(hd, cat, a - hc2_agestart, s) / i_hc.hc_hiv_total(hd, a, s);
            }
          } // end hcTT
        } // end hc1DS
      } // end a
    } // end NS

    for (int s = 0; s < NS; ++s) {
      for (int a = 0; a < p_idx_hiv_first_adult; ++a) {
        for (int hd = 0; hd < hc1DS; ++hd) {
          for (int cat = 0; cat < hcTT; ++cat) {
            if (a < hc2_agestart) {
              auto ltfu_grad = (n_hc.hc1_art_pop(2, hd, a, s) + n_hc.hc1_art_pop(0, hd, a, s)) *
                               p_hc.hc_art_ltfu(t);
              if (i_hc.hc_hiv_total(hd, a, s) > 0) {
                i_hc.art_ltfu_grad(hd, cat, a, s) += ltfu_grad * i_hc.hc_hiv_dist(hd, cat, a, s);
              } else {
                i_hc.art_ltfu_grad(hd, cat, a, s) += ltfu_grad * 0.25;
              }
            } else if (hd < hc2DS) {
              auto ltfu_grad = (n_hc.hc2_art_pop(2, hd, a - hc2_agestart, s) + n_hc.hc2_art_pop(0, hd, a - hc2_agestart, s)) *
                               p_hc.hc_art_ltfu(t);
              if (i_hc.hc_hiv_total(hd, a, s) > 0) {
                i_hc.art_ltfu_grad(hd, cat, a, s) += ltfu_grad * i_hc.hc_hiv_dist(hd, cat, a, s);
              } else {
                i_hc.art_ltfu_grad(hd, cat, a, s) += ltfu_grad * 0.25;
              }
            }
          } // end hcTT
        } // end hc1DS
      } // end a
    } // end NS
  };

  void apply_ltfu_to_hivpop() {
    auto& n_hc = state_next.hc;
    auto& i_hc = intermediate.hc;

    for (int s = 0; s < NS; ++s) {
      for (int a = 0; a < p_idx_hiv_first_adult; ++a) {
        for (int hd = 0; hd < hc1DS; ++hd) {
          for (int cat = 0; cat < hcTT; ++cat) {
            if (a < hc2_agestart) {
              n_hc.hc1_hiv_pop(hd, cat, a, s) += i_hc.art_ltfu_grad(hd, cat, a, s);
            } else if (hd < hc2DS) {
              n_hc.hc2_hiv_pop(hd, cat, a - hc2_agestart, s) += i_hc.art_ltfu_grad(hd, cat, a, s);
            }
          } // end hcTT
        } // end hc1DS
      } // end a
    } // end NS
  };

  void apply_ltfu_to_artpop() {
    auto& n_hc = state_next.hc;
    auto& i_hc = intermediate.hc;

    for (int s = 0; s < NS; ++s) {
      for (int a = 0; a < p_idx_hiv_first_adult; ++a) {
        for (int hd = 0; hd < hc1DS; ++hd) {
          for (int cat = 0; cat < hcTT; ++cat) {
            if (a < hc2_agestart) {
              // hard coded two as this will only occur among children that are on ART more than a year
              n_hc.hc1_art_pop(2, hd, a, s) -= i_hc.art_ltfu_grad(hd, cat, a, s);
            } else if (hd < hc2DS) {
              n_hc.hc2_art_pop(2, hd, a - hc2_agestart, s) -= i_hc.art_ltfu_grad(hd, cat, a, s);
            }
          } // end hcTT
        } // end hc1DS
      } // end a
    } // end NS
  };

  void art_initiation_by_age() {
    const auto& p_hc = pars.hc;
    auto& n_hc = state_next.hc;
    auto& i_hc = intermediate.hc;

    if (p_hc.hc_art_is_age_spec(t)) {
      for (int s = 0; s < NS; ++s) {
        for (int a = 0; a < p_idx_hiv_first_adult; ++a) {
          for (int hd = 0; hd < hc1DS; ++hd) {
            for (int cat = 0; cat < hcTT; ++cat) {
              i_hc.hc_initByAge(hc_age_coarse[a]) += n_hc.hc_art_need_init(hd, cat, a, s) *
                                                          p_hc.hc_art_init_dist(a, t);
            } // end hcTT
          } // end hc1DS
        } // end a
      } // end NS

      for (int ag = 1; ag < hcAG_coarse; ++ag) {
        if (i_hc.hc_initByAge(ag) == 0.0) {
          i_hc.hc_adj(ag) = 1.0;
        } else {
          i_hc.hc_adj(ag) = n_hc.hc_art_init(ag) / i_hc.hc_initByAge(ag);
        }
      }

      for (int s = 0; s < NS; ++s) {
        for (int cat = 0; cat < hcTT; ++cat) {
          for (int a = 0; a < p_idx_hiv_first_adult; ++a) {
            for (int hd = 0; hd < hc1DS; ++hd) {
              auto& coarse_hc_adj = i_hc.hc_adj(hc_age_coarse[a]);
              auto& coarse_hc_art_scalar = i_hc.hc_art_scalar(hc_age_coarse[a]);
              auto hc_art_val_sum = p_hc.hc_art_val(0, t) + p_hc.hc_art_val(1, t) +
                                    p_hc.hc_art_val(2, t) + p_hc.hc_art_val(3, t);
              auto hc_art_val_sum_last = p_hc.hc_art_val(0, t - 1) + p_hc.hc_art_val(1, t - 1) +
                                         p_hc.hc_art_val(2, t - 1) + p_hc.hc_art_val(3, t - 1);
              if (hc_art_val_sum + hc_art_val_sum_last <= 0.0) {
                coarse_hc_art_scalar = 0.0;
              } else {
                coarse_hc_art_scalar = std::min(coarse_hc_adj * p_hc.hc_art_init_dist(a, t), 1.0);
              }

              auto art_initiates = coarse_hc_art_scalar * n_hc.hc_art_need_init(hd, cat, a, s);

              if (a < hc2_agestart) {
                n_hc.hc1_art_pop(0, hd, a, s) += art_initiates;
                n_hc.hc1_hiv_pop(hd, cat, a, s) -= art_initiates;
              } else if (hd < (hc2DS)) {
                n_hc.hc2_art_pop(0, hd, a - hc2_agestart, s) += art_initiates;
                n_hc.hc2_hiv_pop(hd, cat, a - hc2_agestart, s) -=  art_initiates;
              }
            } // end hc1DS
          } // end a
        } // end hcTT
      } // end  NS
    } else {
      for (int s = 0; s < NS; ++s) {
        for (int a = 0; a < p_idx_hiv_first_adult; ++a) {
          for (int hd = 0; hd < hc1DS; ++hd) {
            for (int cat = 0; cat < hcTT; ++cat) {
              i_hc.hc_initByAge(0) += n_hc.hc_art_need_init(hd, cat, a, s) * p_hc.hc_art_init_dist(a, t);
            } // end hcTT
          } // end hc1DS
        } // end a
      } // end  NS

      if (i_hc.hc_initByAge(0) == 0.0) {
        i_hc.hc_adj(0) = 1.0 ;
      } else {
        i_hc.hc_adj(0) =  n_hc.hc_art_init(0) / i_hc.hc_initByAge(0);
      }

      for (int s = 0; s < NS; ++s) {
        for (int cat = 0; cat < hcTT; ++cat) {
          for (int a = 0; a < p_idx_hiv_first_adult; ++a) {
            for (int hd = 0; hd < hc1DS; ++hd) {
              auto hc_art_val_sum = p_hc.hc_art_val(0, t) + p_hc.hc_art_val(0, t - 1);
              if (hc_art_val_sum <= 0) {
                i_hc.hc_art_scalar(0) = 0.0;
              } else {
                i_hc.hc_art_scalar(0) = std::min(i_hc.hc_adj(0) * p_hc.hc_art_init_dist(a, t), 1.0);
              }

              auto art_initiates = i_hc.hc_art_scalar(0) * n_hc.hc_art_need_init(hd, cat, a, s);
              if (a < hc2_agestart) {
                n_hc.hc1_art_pop(0, hd, a, s) += art_initiates;
                n_hc.hc1_hiv_pop(hd, cat, a, s) -= art_initiates;
              } else if (hd < (hc2DS)) {
                n_hc.hc2_art_pop(0, hd, a - hc2_agestart, s) += art_initiates;
                n_hc.hc2_hiv_pop(hd, cat, a - hc2_agestart, s) -= art_initiates;
              }
            } // end hc1DS
          } // end a
        } // end hcTT
      } // end  NS
    } // end if
  };

  void fill_total_pop_outputs() {
    auto& n_ha = state_next.ha;
    auto& n_hc = state_next.hc;

    for (int hd = 0; hd < hDS; ++hd) {
      for (int a = 0; a < p_idx_hiv_first_adult; ++a) {
        for (int s = 0; s < NS; ++s) {
          for (int cat = 0; cat < hcTT; ++cat) {
            if (a < hc2_agestart) {
              n_ha.p_hiv_deaths(a,s) += n_hc.hc1_noart_aids_deaths(hd, cat, a, s);
            } else if (hd < hc2DS) {
              n_ha.p_hiv_deaths(a,s) += n_hc.hc2_noart_aids_deaths(hd, cat, a - hc2_agestart, s);
            }
          } // end hcTT

          for (int dur = 0; dur < hTS; ++dur) {
            if (a < hc2_agestart) {
              n_ha.p_hiv_deaths(a,s) += n_hc.hc1_art_aids_deaths(dur, hd, a, s);
            } else if (hd < hc2DS) {
              n_ha.p_hiv_deaths(a,s) += n_hc.hc2_art_aids_deaths(dur, hd, a - hc2_agestart, s);
            }
          } // end dur
        } // end NS
      } // end a
    } // end hDS

    for (int a = 0; a < p_idx_hiv_first_adult; ++a) {
      for (int s = 0; s < NS; ++s) {
        n_ha.p_hiv_pop(a, s) += n_ha.p_infections(a, s);
        n_ha.p_hiv_pop(a, s) -= n_ha.p_hiv_deaths(a, s);
      }
    }
  };
};

}
}
