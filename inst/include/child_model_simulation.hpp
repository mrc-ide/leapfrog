#pragma once

#include "intermediate_data.hpp"

namespace leapfrog {

namespace internal {

template<typename ModelVariant, typename real_type>
void run_child_ageing(int t,
                      const Parameters<ModelVariant, real_type> &pars,
                      const State<ModelVariant, real_type> &state_curr,
                      State<ModelVariant, real_type> &state_next,
                      IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "run_hiv_child_infections can only be called for model variants where run_child_model is true");
  constexpr auto ss_d = StateSpace<ModelVariant>().dp;
  constexpr auto ss_h = StateSpace<ModelVariant>().hiv;
  constexpr auto ss_c = StateSpace<ModelVariant>().children;
  const auto& p_dm = pars.dp.demography;
  const auto& p_hc = pars.children.children;
  const auto& p_op = pars.options;
  const auto& c_hc = state_curr.children;
  auto& n_hc = state_next.children;

  for (int s = 0; s < ss_d.NS; ++s) {
    // less than 5 because there is a cd4 transition between ages 4 and 5
    for (int a = 1; a < ss_c.hc2_agestart; ++a) {
      for (int hd = 0; hd < ss_c.hc1DS; ++hd) {
        for (int cat = 0; cat < ss_c.hcTT; ++cat) {
          n_hc.hc1_hiv_pop(hd, cat, a, s) += c_hc.hc1_hiv_pop(hd, cat, a - 1, s) * p_dm.survival_probability(a, s, t);
        }
        for (int dur = 0; dur < ss_h.hTS; ++dur) {
          n_hc.hc1_art_pop(dur, hd, a, s) += c_hc.hc1_art_pop(dur, hd, a - 1, s) * p_dm.survival_probability(a, s, t);
        }
      }
    }
  }

  for (int s = 0; s < ss_d.NS; ++s) {
    for (int hd = 0; hd < ss_c.hc1DS; ++hd) {
      for (int hd_alt = 0; hd_alt < ss_c.hc2DS; ++hd_alt) {
        for (int cat = 0; cat < ss_c.hcTT; ++cat) {
          n_hc.hc2_hiv_pop(hd_alt, cat, 0, s) += c_hc.hc1_hiv_pop(hd, cat, ss_c.hc1_ageend, s) *
                                                 p_dm.survival_probability(ss_c.hc2_agestart, s, t) *
                                                 p_hc.hc_cd4_transition(hd_alt, hd);
        }
        for (int dur = 0; dur < ss_h.hTS; ++dur) {
          n_hc.hc2_art_pop(dur, hd_alt, 0, s) += c_hc.hc1_art_pop(dur, hd, ss_c.hc1_ageend, s) *
                                                 p_dm.survival_probability(ss_c.hc2_agestart, s, t) *
                                                 p_hc.hc_cd4_transition(hd_alt, hd);
        }
      }
    }
  }

  for (int s = 0; s < ss_d.NS; ++s) {
    for (int a = (ss_c.hc2_agestart + 1); a < p_op.p_idx_fertility_first; ++a) {
      for (int hd = 0; hd < ss_c.hc2DS; ++hd) {
        for (int cat = 0; cat < ss_c.hcTT; ++cat) {
          n_hc.hc2_hiv_pop(hd, cat, a - ss_c.hc2_agestart, s) += c_hc.hc2_hiv_pop(hd, cat, a - ss_c.hc2_agestart - 1, s) *
                                                                 p_dm.survival_probability(a, s, t);
        }
        for (int dur = 0; dur < ss_h.hTS; ++dur) {
          n_hc.hc2_art_pop(dur, hd, a - ss_c.hc2_agestart, s) += c_hc.hc2_art_pop(dur, hd, a - ss_c.hc2_agestart - 1, s) *
                                                                 p_dm.survival_probability(a, s, t);
        }
      }
    }
  }

}

template<typename ModelVariant, typename real_type>
void run_wlhiv_births(int t,
                      const Parameters<ModelVariant, real_type> &pars,
                      const State<ModelVariant, real_type> &state_curr,
                      State<ModelVariant, real_type> &state_next,
                      IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "run_wlhiv_births can only be called for model variants where run_child_model is true");
  constexpr auto ss_h = StateSpace<ModelVariant>().hiv;
  const auto& p_hc = pars.children.children;
  const auto& p_dm = pars.dp.demography;
  const auto& p_op = pars.options;
  const auto& c_ba = state_curr.hiv;
  auto& n_ha = state_next.hiv;
  auto& n_dp = state_next.dp;
  auto& n_hc = state_next.children;
  auto& i_hc = intermediate.children;

  i_hc.asfr_sum = 0.0;
  for (int a = 0; a < p_op.p_fertility_age_groups; ++a) {
    i_hc.asfr_sum += p_dm.age_specific_fertility_rate(a, t);
  } // end a

  for (int a = 0; a < p_op.p_fertility_age_groups; ++a) {
    i_hc.nHIVcurr = 0.0;
    i_hc.nHIVlast = 0.0;
    i_hc.df = 0.0;

    for (int hd = 0; hd < ss_h.hDS; ++hd) {
      i_hc.nHIVcurr += n_ha.h_hiv_adult(hd, a, 1);
      i_hc.nHIVlast += c_ba.h_hiv_adult(hd, a, 1);
      for (int ht = 0; ht < ss_h.hTS; ++ht) {
        i_hc.nHIVcurr += n_ha.h_art_adult(ht, hd, a, 1);
        i_hc.nHIVlast += c_ba.h_art_adult(ht, hd, a, 1);
      } // end hTS
    } // end hDS

    i_hc.prev = i_hc.nHIVcurr / n_dp.p_total_pop(a + 15, 1);

    for (int hd = 0; hd < ss_h.hDS; ++hd) {
      i_hc.df += p_hc.local_adj_factor * p_hc.fert_mult_by_age(a) * p_hc.fert_mult_off_art(hd) *
                 (n_ha.h_hiv_adult(hd, a, 1) + c_ba.h_hiv_adult(hd, a, 1)) / 2;
      // women on ART less than 6 months use the off art fertility multiplier
      i_hc.df += p_hc.local_adj_factor * p_hc.fert_mult_by_age(a) * p_hc.fert_mult_off_art(hd) *
                 (n_ha.h_art_adult(0, hd, a, 1) + c_ba.h_art_adult(0, hd, a, 1)) / 2;
      for (int ht = 1; ht < ss_h.hTS; ++ht) {
        i_hc.df += p_hc.local_adj_factor * p_hc.fert_mult_on_art(a) *
                   (n_ha.h_art_adult(ht, hd, a, 1) + c_ba.h_art_adult(ht, hd, a, 1)) / 2;
      } // end hTS
    } // end hDS


    if (i_hc.nHIVcurr > 0) {
      auto midyear_fertileHIV = (i_hc.nHIVcurr + i_hc.nHIVlast) / 2;
      i_hc.df = i_hc.df / midyear_fertileHIV;
      i_hc.birthsCurrAge = midyear_fertileHIV * p_hc.total_fertility_rate(t) *
        i_hc.df / (i_hc.df * i_hc.prev + 1 - i_hc.prev) *
        p_dm.age_specific_fertility_rate(a, t) / i_hc.asfr_sum ;
    } else {
      auto midyear_fertileHIV = (i_hc.nHIVcurr + i_hc.nHIVlast) / 2;
      i_hc.df = 1;
      i_hc.birthsCurrAge = midyear_fertileHIV * p_hc.total_fertility_rate(t) *
        i_hc.df / (i_hc.df * i_hc.prev + 1 - i_hc.prev) *
        p_dm.age_specific_fertility_rate(a, t) / i_hc.asfr_sum ;
    }

    i_hc.birthsHE += i_hc.birthsCurrAge;
    if (a < 9) {
      i_hc.births_HE_15_24 += i_hc.birthsCurrAge;
    }
  } // end a
  n_hc.hiv_births = i_hc.birthsHE;
}

template<typename ModelVariant, typename real_type>
void run_wlhiv_births_input_mat_prev(int t,
                                     const Parameters<ModelVariant, real_type> &pars,
                                     const State<ModelVariant, real_type> &state_curr,
                                     State<ModelVariant, real_type> &state_next,
                                     IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "run_wlhiv_births_input_mat_prev can only be called for model variants where run_child_model is true");
  const auto& p_hc = pars.children.children;
  auto& n_hc = state_next.children;

  n_hc.hiv_births = p_hc.mat_hiv_births(t);
}

template<typename ModelVariant, typename real_type>
void calc_hiv_negative_pop(int t,
                           const Parameters<ModelVariant, real_type> &pars,
                           const State<ModelVariant, real_type> &state_curr,
                           State<ModelVariant, real_type> &state_next,
                           IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "calc_hiv_negative_pop can only be called for model variants where run_child_model is true");
  constexpr auto ss_d = StateSpace<ModelVariant>().dp;
  constexpr auto ss_h = StateSpace<ModelVariant>().hiv;
  const auto& p_dm = pars.dp.demography;
  auto& n_ha = state_next.hiv;
  auto& n_dp = state_next.dp;
  auto& i_hc = intermediate.children;

  for (int s = 0; s < ss_d.NS; ++s) {
    for (int a = 0; a < ss_h.hAG; ++a) {
      i_hc.p_hiv_neg_pop(a, s) = n_dp.p_total_pop(a, s) - n_ha.p_hiv_pop(a, s);
    }// end a
  }// end s
}

template<typename ModelVariant, typename real_type>
void adjust_hiv_births(int t,
                       const Parameters<ModelVariant, real_type> &pars,
                       const State<ModelVariant, real_type> &state_curr,
                       State<ModelVariant, real_type> &state_next,
                       IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "adjust_hiv_births can only be called for model variants where run_child_model is true");
  const auto& p_hc = pars.children.children;
  auto& n_hc = state_next.children;

  if (p_hc.abortion(t, 1) == 1) {
    n_hc.hiv_births -= n_hc.hiv_births * p_hc.abortion(t, 0);
  } else {
    n_hc.hiv_births -=  p_hc.abortion(t, 0);
  }
}

template<typename ModelVariant, typename real_type>
void convert_PMTCT_num_to_perc(int t,
                               const Parameters<ModelVariant, real_type> &pars,
                               const State<ModelVariant, real_type> &state_curr,
                               State<ModelVariant, real_type> &state_next,
                               IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "convert_PMTCT_num_to_perc can only be called for model variants where run_child_model is true");
  constexpr auto ss_c = StateSpace<ModelVariant>().children;
  const auto& p_hc = pars.children.children;
  auto& n_hc = state_next.children;
  auto& i_hc = intermediate.children;

  i_hc.sumARV = 0.0;
  for (int hp = 0; hp < ss_c.hPS; ++hp) {
    i_hc.sumARV += p_hc.PMTCT(hp, t);
  }

  i_hc.need_PMTCT = std::max(i_hc.sumARV, n_hc.hiv_births);

  i_hc.OnPMTCT = i_hc.sumARV + p_hc.patients_reallocated(t);
  i_hc.OnPMTCT = std::min(i_hc.OnPMTCT, i_hc.need_PMTCT);

  // replace all instances of coverage input as numbers with percentage covered
  if (p_hc.PMTCT_input_is_percent(t)) {
    for (int hp = 0; hp < ss_c.hPS; ++hp) {
      i_hc.PMTCT_coverage(hp) = p_hc.PMTCT(hp, t) / 100;
    } // end hPS
    i_hc.sumARV = i_hc.sumARV * i_hc.need_PMTCT;
  } else {
    for (int hp = 0; hp < ss_c.hPS; ++hp) {
      if (i_hc.sumARV == 0) {
        i_hc.PMTCT_coverage(hp) = 0.0;
      } else {
        i_hc.PMTCT_coverage(hp) = p_hc.PMTCT(hp, t) / i_hc.sumARV *
                                  i_hc.OnPMTCT / i_hc.need_PMTCT;
      }
    } // end hPS
  }

  i_hc.PMTCT_coverage(4) = i_hc.PMTCT_coverage(4) * p_hc.PMTCT_dropout(0, t);
  i_hc.PMTCT_coverage(5) = i_hc.PMTCT_coverage(5) * p_hc.PMTCT_dropout(1, t);
}

template<typename ModelVariant, typename real_type>
void convert_PMTCT_pre_bf(int t,
                          const Parameters<ModelVariant, real_type> &pars,
                          const State<ModelVariant, real_type> &state_curr,
                          State<ModelVariant, real_type> &state_next,
                          IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "convert_PMTCT_num_to_perc can only be called for model variants where run_child_model is true");
  constexpr auto ss_c = StateSpace<ModelVariant>().children;
  const auto& p_hc = pars.children.children;
  auto& i_hc = intermediate.children;

  // TODO: Maggie to confirm why Option A/B alt tr aren't used
  for (int hp = 0; hp < ss_c.hPS; ++hp) {
    i_hc.PMTCT_coverage(hp) *= 1.0 - p_hc.PMTCT_transmission_rate(4, hp, 0);
  } // end hPS
}

template<typename ModelVariant, typename real_type>
void calc_wlhiv_cd4_proportion(int t,
                               const Parameters<ModelVariant, real_type> &pars,
                               const State<ModelVariant, real_type> &state_curr,
                               State<ModelVariant, real_type> &state_next,
                               IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "calc_wlhiv_cd4_proportion can only be called for model variants where run_child_model is true");
  const auto& p_hc = pars.children.children;
  auto& n_ha = state_next.hiv;
  auto& i_hc = intermediate.children;

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
}

template<typename ModelVariant, typename real_type>
void adjust_option_A_B_tr(int t,
                          const Parameters<ModelVariant, real_type> &pars,
                          const State<ModelVariant, real_type> &state_curr,
                          State<ModelVariant, real_type> &state_next,
                          IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "adjust_option_A_B_tr can only be called for model variants where run_child_model is true");
  const auto& p_hc = pars.children.children;
  auto& i_hc = intermediate.children;

  // Option A and B were only authorized for women with greater than 350 CD4, so if the percentage of women
  // on option A/B > the proportion of women in this cd4 category, we assume that some must have a cd4 less than 350
  // option AB will be less effective for these women so we adjust for that
  internal::calc_wlhiv_cd4_proportion(t, pars, state_curr, state_next, intermediate);

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
}

template<typename ModelVariant, typename real_type>
void adjust_option_A_B_bf_tr(int t,
                             const Parameters<ModelVariant, real_type> &pars,
                             const State<ModelVariant, real_type> &state_curr,
                             State<ModelVariant, real_type> &state_next,
                             IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "adjust_option_A_B_bf_tr can only be called for model variants where run_child_model is true");
  const auto& p_hc = pars.children.children;
  auto& i_hc = intermediate.children;

  internal::calc_wlhiv_cd4_proportion(t, pars, state_curr, state_next, intermediate);

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
}

template<typename ModelVariant, typename real_type>
void maternal_incidence_in_pregnancy_tr(int t,
                                        const Parameters<ModelVariant, real_type> &pars,
                                        const State<ModelVariant, real_type> &state_curr,
                                        State<ModelVariant, real_type> &state_next,
                                        IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "maternal_incidence_in_pregnancy_tr can only be called for model variants where run_child_model is true");
  const auto& p_dm = pars.dp.demography;
  const auto& p_op = pars.options;
  const auto& p_hc = pars.children.children;
  auto& n_dp = state_next.dp;
  auto& n_ha = state_next.hiv;
  auto& i_hc = intermediate.children;

  // Transmission due to incident infections
  i_hc.asfr_sum = 0.0;
  for (int a = 0; a < p_op.p_fertility_age_groups; ++a) {
    i_hc.asfr_sum += p_dm.age_specific_fertility_rate(a, t);
  } // end a

  if (p_hc.mat_prev_input(t)) {
    for (int a = 0; a < p_op.p_fertility_age_groups; ++a) {
      auto asfr_weight = p_dm.age_specific_fertility_rate(a, t) / i_hc.asfr_sum;
      i_hc.age_weighted_hivneg += asfr_weight * p_hc.adult_female_hivnpop(a, t); // HIV negative 15-49 women weighted for ASFR
      i_hc.age_weighted_infections += asfr_weight * p_hc.adult_female_infections(a, t); // newly infected 15-49 women, weighted for ASFR
    } // end a

    if (i_hc.age_weighted_hivneg > 0.0) {
      i_hc.incidence_rate_wlhiv = i_hc.age_weighted_infections / i_hc.age_weighted_hivneg;
      // 0.75 is 9/12, gestational period, index 7 in the vertical trasnmission object is the index for maternal seroconversion
      i_hc.perinatal_transmission_from_incidence = i_hc.incidence_rate_wlhiv * 0.75 *
                                                   (p_hc.total_births(t) - i_hc.need_PMTCT) *
                                                   p_hc.vertical_transmission_rate(7, 0);
    } else {
      i_hc.incidence_rate_wlhiv = 0.0;
      i_hc.perinatal_transmission_from_incidence = 0.0;
    }
  } else {
    for (int a = 0; a < p_op.p_fertility_age_groups; ++a) {
      auto asfr_weight = p_dm.age_specific_fertility_rate(a, t) / i_hc.asfr_sum;
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
}

template<typename ModelVariant, typename real_type>
void perinatal_tr(int t,
                  const Parameters<ModelVariant, real_type> &pars,
                  const State<ModelVariant, real_type> &state_curr,
                  State<ModelVariant, real_type> &state_next,
                  IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "perinatal_tr can only be called for model variants where run_child_model is true");
  const auto& p_hc = pars.children.children;
  auto& n_dp = state_next.dp;
  auto& i_hc = intermediate.children;

  i_hc.births_sum = n_dp.births;

  // TODO: add in patients reallocated
  internal::convert_PMTCT_num_to_perc(t, pars, state_curr, state_next, intermediate);
  internal::adjust_option_A_B_tr(t, pars, state_curr, state_next, intermediate);
  internal::calc_hiv_negative_pop(t, pars, state_curr, state_next, intermediate);

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
  internal::maternal_incidence_in_pregnancy_tr(t, pars, state_curr, state_next, intermediate);
  if (i_hc.need_PMTCT > 0.0) {
    i_hc.perinatal_transmission_rate += i_hc.perinatal_transmission_from_incidence / i_hc.need_PMTCT;
  }

}

template<typename ModelVariant, typename real_type>
void maternal_incidence_in_bf_tr(int t,
                                 const Parameters<ModelVariant, real_type> &pars,
                                 const State<ModelVariant, real_type> &state_curr,
                                 State<ModelVariant, real_type> &state_next,
                                 IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "maternal_incidence_in_bf_tr can only be called for model variants where run_child_model is true");
  constexpr auto ss_c = StateSpace<ModelVariant>().children;
  const auto& p_hc = pars.children.children;
  auto& i_hc = intermediate.children;

  for (int bf = 0; bf < ss_c.hBF; ++bf) {
    i_hc.bf_at_risk += i_hc.incidence_rate_wlhiv / 12 * 2 *
                       (1 - p_hc.breastfeeding_duration_no_art(bf, t));
  }
   i_hc.bf_incident_hiv_transmission_rate = i_hc.bf_at_risk * p_hc.vertical_transmission_rate(7, 1);
}

template<typename ModelVariant, typename real_type>
void bf_dropout(int t,
                const Parameters<ModelVariant, real_type> &pars,
                const State<ModelVariant, real_type> &state_curr,
                State<ModelVariant, real_type> &state_next,
                IntermediateData<ModelVariant, real_type> &intermediate,
                int bf) {
  static_assert(ModelVariant::run_child_model,
                "bf_dropout can only be called for model variants where run_child_model is true");
  constexpr auto ss_c = StateSpace<ModelVariant>().children;
  const auto& p_hc = pars.children.children;
  auto& i_hc = intermediate.children;

  if(bf < 6){
    i_hc.PMTCT_coverage(0) *= (1 - p_hc.PMTCT_dropout(4, t) * 2); //opt A
    i_hc.PMTCT_coverage(1) *= (1 - p_hc.PMTCT_dropout(4, t) * 2); //opt B
    i_hc.PMTCT_coverage(4) *= (1 - p_hc.PMTCT_dropout(4, t) * 2); //before pregnancy
    i_hc.PMTCT_coverage(5) *= (1 - p_hc.PMTCT_dropout(4, t) * 2); //>4 weeks
    i_hc.PMTCT_coverage(6) *= (1 - p_hc.PMTCT_dropout(4, t) * 2); //<4 weeks
  }else{
    i_hc.PMTCT_coverage(0) *= (1 - p_hc.PMTCT_dropout(5, t) * 2); //opt A
    i_hc.PMTCT_coverage(1) *= (1 - p_hc.PMTCT_dropout(5, t) * 2); //opt B
    i_hc.PMTCT_coverage(4) *= (1 - p_hc.PMTCT_dropout(5, t) * 2); //before pregnancy
    i_hc.PMTCT_coverage(5) *= (1 - p_hc.PMTCT_dropout(5, t) * 2); //>4 weeks
    i_hc.PMTCT_coverage(6) *= (1 - p_hc.PMTCT_dropout(5, t) * 2); //<4 weeks
  }

}

template<typename ModelVariant, typename real_type>
void run_bf_transmission_rate(int t,
                              const Parameters<ModelVariant, real_type> &pars,
                              const State<ModelVariant, real_type> &state_curr,
                              State<ModelVariant, real_type> &state_next,
                              IntermediateData<ModelVariant, real_type> &intermediate,
                              int bf_start, int bf_end, int index) {
  static_assert(ModelVariant::run_child_model,
                "run_bf_transmission_rate can only be called for model variants where run_child_model is true");
  constexpr auto ss_c = StateSpace<ModelVariant>().children;
  const auto& p_hc = pars.children.children;
  auto& i_hc = intermediate.children;
  auto& n_hc = state_next.children;

  for (int bf = bf_start; bf < bf_end; bf++) {
    if(bf == 0){
      //Perinatal transmission accounts for transmission up to 6 weeks, so we only use 1/4 of
      //transmission from the first breastfeeding period
      i_hc.bf_scalar = 0.25;
    }else{
      i_hc.bf_scalar = 1.0;
      //dropout only occurs after the first month of breastfeeding
      internal::bf_dropout(t, pars, state_curr, state_next, intermediate, bf);
    }
    // i_hc.perinatal_transmission_rate_bf_calc is the transmission that has already occurred due to perinatal transmission
    // i_hc.percent_no_treatment is the percentage of women who are still vulnerable to HIV transmission to their babies
    i_hc.percent_no_treatment = 1 - i_hc.perinatal_transmission_rate_bf_calc - i_hc.bf_transmission_rate(index) ;
    if(index > 0){
      for (int bf = 0; bf < index; ++bf) {
        i_hc.percent_no_treatment -= i_hc.bf_transmission_rate(bf);
      }
    }

    for (int hp = 0; hp < ss_c.hPS; hp++) {
      i_hc.percent_on_treatment = 0;
      i_hc.percent_no_treatment -=  i_hc.PMTCT_coverage(hp);

      if(hp > 1){
        if(hp == 2){
          auto tr = i_hc.PMTCT_coverage(hp) *
            //sdnvp stratifies transmission by CD4, but spectrum only uses one
            p_hc.PMTCT_transmission_rate(0, hp, 1) *
            2 * (1 - p_hc.breastfeeding_duration_art(bf, t)) * i_hc.bf_scalar;
          i_hc.PMTCT_coverage(hp) -= tr;
          i_hc.bf_transmission_rate(index) += tr;
        }else{
          auto tr = i_hc.PMTCT_coverage(hp) *
            p_hc.PMTCT_transmission_rate(4, hp, 1) *
            2 * (1 - p_hc.breastfeeding_duration_art(bf, t)) * i_hc.bf_scalar;
          i_hc.PMTCT_coverage(hp) -= tr;
          i_hc.bf_transmission_rate(index) += tr;
        }
      }
    }

    // No treatment
    if (p_hc.breastfeeding_duration_no_art(bf, t) < 1) {
      i_hc.percent_no_treatment = std::max(i_hc.percent_no_treatment, 0.0);
      auto untreated_vertical_bf_tr = i_hc.prop_wlhiv_lt200 * p_hc.vertical_transmission_rate(4, 1) +
        i_hc.prop_wlhiv_200to350 * p_hc.vertical_transmission_rate(2, 1) +
        i_hc.prop_wlhiv_gte350 * p_hc.vertical_transmission_rate(0, 1);
      i_hc.bf_transmission_rate(index) += i_hc.bf_scalar *
        i_hc.percent_no_treatment *
        untreated_vertical_bf_tr *
        2 * (1 - p_hc.breastfeeding_duration_no_art(bf, t));
    }
  }

}

template<typename ModelVariant, typename real_type>
void nosocomial_infections(int t,
                           const Parameters<ModelVariant, real_type> &pars,
                           const State<ModelVariant, real_type> &state_curr,
                           State<ModelVariant, real_type> &state_next,
                           IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "nosocomial_infections can only be called for model variants where run_child_model is true");
  constexpr auto ss_d = StateSpace<ModelVariant>().dp;
  constexpr auto ss_c = StateSpace<ModelVariant>().children;
  const auto& p_hc = pars.children.children;
  auto& n_ha = state_next.hiv;
  auto& n_hc = state_next.children;

  for (int s = 0; s < ss_d.NS; ++s) {
    // Run only first 5 age groups in total population 0, 1, 2, 3, 4
    for (int a = 0; a < ss_c.hc2_agestart; ++a) {
      if (p_hc.hc_nosocomial(t) > 0) {
        // 5.0 is used because we want to evenly distribute across the 5 age groups in 0-4
        n_ha.p_infections(a, s) = p_hc.hc_nosocomial(t) / (5.0 * ss_d.NS);
        // Putting all nosocomial acquired HIV infections in perinatally acquired infection timing and highest CD4 category to match Spectrum implementation
        n_hc.hc1_hiv_pop(0, 0, a, s) += n_ha.p_infections(a, s);
      }
    } // end a
  } // end NS
}

template<typename ModelVariant, typename real_type>
void add_infections(int t,
                    const Parameters<ModelVariant, real_type> &pars,
                    const State<ModelVariant, real_type> &state_curr,
                    State<ModelVariant, real_type> &state_next,
                    IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "add_infections can only be called for model variants where run_child_model is true");
  constexpr auto ss_d = StateSpace<ModelVariant>().dp;
  constexpr auto ss_c = StateSpace<ModelVariant>().children;
  const auto& p_dm = pars.dp.demography;
  const auto& p_hc = pars.children.children;
  auto& n_ha = state_next.hiv;
  auto& n_dp = state_next.dp;
  auto& n_hc = state_next.children;
  auto& i_hc = intermediate.children;

  if (n_hc.hiv_births > 0) {
    internal::perinatal_tr(t, pars, state_curr, state_next, intermediate);

    // Perinatal transmission
    auto perinatal_transmission_births = n_hc.hiv_births * i_hc.perinatal_transmission_rate;
    for (int s = 0; s < ss_d.NS; ++s) {
      for (int hd = 0; hd < ss_c.hc1DS; ++hd) {
      n_hc.hc1_hiv_pop(hd, 0, 0, s) += perinatal_transmission_births *
                                           p_dm.births_sex_prop(s, t) * p_hc.hc1_cd4_dist(hd);
      } // end hc1DS
      auto perinatal_births_by_sex = perinatal_transmission_births * p_dm.births_sex_prop(s, t);
      n_ha.p_infections(0, s) += perinatal_births_by_sex;
      n_hc.infection_by_type(0,0,s) += perinatal_births_by_sex;
    } // end NS

    // Breastfeeding transmission
    // 0-6
    internal::maternal_incidence_in_bf_tr(t, pars, state_curr, state_next, intermediate);
    internal::adjust_option_A_B_bf_tr(t, pars, state_curr, state_next, intermediate);
    internal::convert_PMTCT_pre_bf(t, pars, state_curr, state_next, intermediate);
    internal::run_bf_transmission_rate(t, pars, state_curr, state_next, intermediate, 0, 3, 0);
    n_hc.hc_stacked_bar(6,1) += i_hc.bf_incident_hiv_transmission_rate;
    real_type total_births = 0.0;
    if (p_hc.mat_prev_input(t)) {
      total_births = p_hc.total_births(t);
    } else {
      total_births = n_dp.births;
    }

    // 0-6
    for (int s = 0; s < ss_d.NS; ++s) {
      auto bf_hiv_by_sex = n_hc.hiv_births * p_dm.births_sex_prop(s, t) *
                           i_hc.bf_transmission_rate(0);
      // vertical infection from maternal infection during breastfeeding
      bf_hiv_by_sex += (total_births - n_hc.hiv_births) * p_dm.births_sex_prop(s, t) *
                       i_hc.bf_incident_hiv_transmission_rate;
      for (int hd = 0; hd < ss_c.hc1DS; ++hd) {
        n_hc.hc1_hiv_pop(hd, 1, 0, s) += p_hc.hc1_cd4_dist(hd) * bf_hiv_by_sex;
      } // end hc1DS
      n_ha.p_infections(0, s) += bf_hiv_by_sex;
      n_hc.infection_by_type(1,0,s) += bf_hiv_by_sex;
    } // end NS

    // 6-12
    internal::run_bf_transmission_rate(t, pars, state_curr, state_next, intermediate, 3, 6, 1);
    for (int s = 0; s < ss_d.NS; ++s) {
      auto bf_hiv_by_sex = n_hc.hiv_births * p_dm.births_sex_prop(s, t) * i_hc.bf_transmission_rate(1);
      for (int hd = 0; hd < ss_c.hc1DS; ++hd) {
        n_hc.hc1_hiv_pop(hd, 2, 0, s) += p_hc.hc1_cd4_dist(hd) * bf_hiv_by_sex;
      } // end hc1DS
      n_ha.p_infections(0, s) += bf_hiv_by_sex;
      n_hc.infection_by_type(2,0,s) += bf_hiv_by_sex;
    } // end NS

    // 12 plus
    internal::run_bf_transmission_rate(t, pars, state_curr, state_next, intermediate, 6, 12, 2);
    internal::run_bf_transmission_rate(t, pars, state_curr, state_next, intermediate, 12, ss_c.hBF, 3);
    for (int s = 0; s < ss_d.NS; ++s) {
      for (int hd = 0; hd < ss_c.hc1DS; ++hd) {
        if(s == 0){
          auto uninfected_prop_12_24 = (n_dp.p_total_pop(1, s) - n_ha.p_hiv_pop(1, s)) /
            (n_dp.p_total_pop(1, 0) - n_ha.p_hiv_pop(1, 0) + n_dp.p_total_pop(1, 1) - n_ha.p_hiv_pop(1, 1));

          auto uninfected_prop_24_plus = (n_dp.p_total_pop(2, s) - n_ha.p_hiv_pop(2, s)) /
            (n_dp.p_total_pop(2, 0) - n_ha.p_hiv_pop(2, 0) + n_dp.p_total_pop(2, 1) - n_ha.p_hiv_pop(2, 1));


          auto bf_hiv_transmission_12_24 = n_hc.hiv_births * i_hc.bf_transmission_rate(2) * uninfected_prop_12_24;
          auto bf_hiv_transmission_24_plus = n_hc.hiv_births * i_hc.bf_transmission_rate(3) * uninfected_prop_24_plus;
          // 12-24
          n_hc.hc1_hiv_pop(hd, 3, 1, s) += p_hc.hc1_cd4_dist(hd) * bf_hiv_transmission_12_24;
          n_ha.p_infections(1, s) += p_hc.hc1_cd4_dist(hd) * bf_hiv_transmission_12_24;
          n_hc.infection_by_type(3,1,s) += p_hc.hc1_cd4_dist(hd) * bf_hiv_transmission_12_24;

          // 24 plus
          n_hc.hc1_hiv_pop(hd, 3, 2, s) += p_hc.hc1_cd4_dist(hd) * bf_hiv_transmission_24_plus;
          n_ha.p_infections(2, s) += p_hc.hc1_cd4_dist(hd) * bf_hiv_transmission_24_plus;
          n_hc.infection_by_type(3,2,s) += p_hc.hc1_cd4_dist(hd) * bf_hiv_transmission_24_plus;
        }else{
          auto uninfected_prop_12_24 = (n_dp.p_total_pop(1, s) - n_ha.p_hiv_pop(1, s)) /
            (n_dp.p_total_pop(1, 0) - n_ha.p_hiv_pop(1, 0) - n_ha.p_infections(1, 0) + n_dp.p_total_pop(1, 1) - n_ha.p_hiv_pop(1, 1));

          auto uninfected_prop_24_plus = (n_dp.p_total_pop(2, s) - n_ha.p_hiv_pop(2, s)) /
            (n_dp.p_total_pop(2, 0) - n_ha.p_hiv_pop(2, 0) - n_ha.p_infections(2, 0) + n_dp.p_total_pop(2, 1) - n_ha.p_hiv_pop(2, 1));


          auto bf_hiv_transmission_12_24 = n_hc.hiv_births * i_hc.bf_transmission_rate(2) * uninfected_prop_12_24;
          auto bf_hiv_transmission_24_plus = n_hc.hiv_births * i_hc.bf_transmission_rate(3) * uninfected_prop_24_plus;
          // 12-24
          n_hc.hc1_hiv_pop(hd, 3, 1, s) += p_hc.hc1_cd4_dist(hd) * bf_hiv_transmission_12_24;
          n_ha.p_infections(1, s) += p_hc.hc1_cd4_dist(hd) * bf_hiv_transmission_12_24;
          n_hc.infection_by_type(3,1,s) += p_hc.hc1_cd4_dist(hd) * bf_hiv_transmission_12_24;

          // 24 plus
          n_hc.hc1_hiv_pop(hd, 3, 2, s) += p_hc.hc1_cd4_dist(hd) * bf_hiv_transmission_24_plus;
          n_ha.p_infections(2, s) += p_hc.hc1_cd4_dist(hd) * bf_hiv_transmission_24_plus;
          n_hc.infection_by_type(3,2,s) += p_hc.hc1_cd4_dist(hd) * bf_hiv_transmission_24_plus;
        }
      } // end hc1DS
    } // end NS
  }
}

template<typename ModelVariant, typename real_type>
void art_eligibility_by_age(int t,
                            const Parameters<ModelVariant, real_type> &pars,
                            const State<ModelVariant, real_type> &state_curr,
                            State<ModelVariant, real_type> &state_next,
                            IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "art_eligibility_by_age can only be called for model variants where run_child_model is true");
  constexpr auto ss_d = StateSpace<ModelVariant>().dp;
  constexpr auto ss_c = StateSpace<ModelVariant>().children;
  const auto& p_hc = pars.children.children;
  auto& n_hc = state_next.children;

  // all children under a certain age eligible for ART
  for (int s = 0; s < ss_d.NS; ++s) {
    for (int a = 0; a < p_hc.hc_art_elig_age(t); ++a) {
      for (int cat = 0; cat < ss_c.hcTT; ++cat) {
        for (int hd = 0; hd < ss_c.hc1DS; ++hd) {
          if (a < ss_c.hc2_agestart) {
            n_hc.hc_art_need_init(hd, cat, a, s) += n_hc.hc1_hiv_pop(hd, cat, a, s);
          } else if (hd < ss_c.hc2DS) {
            n_hc.hc_art_need_init(hd, cat, a, s) += n_hc.hc2_hiv_pop(hd, cat, a - ss_c.hc2_agestart, s);
          }
        } // end ss_c.hc1DS
      } // end a
    } // end hcTT
  } // end ss_d.NS
}

template<typename ModelVariant, typename real_type>
void art_eligibility_by_cd4(int t,
                            const Parameters<ModelVariant, real_type> &pars,
                            const State<ModelVariant, real_type> &state_curr,
                            State<ModelVariant, real_type> &state_next,
                            IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "art_eligibility_by_cd4 can only be called for model variants where run_child_model is true");
  constexpr auto ss_d = StateSpace<ModelVariant>().dp;
  constexpr auto ss_c = StateSpace<ModelVariant>().children;
  const auto& p_hc = pars.children.children;
  const auto& p_op = pars.options;

  // all children under a certain CD4 eligible for ART
  for (int s = 0; s < ss_d.NS; ++s) {
    for (int cat = 0; cat < ss_c.hcTT; ++cat) {
      for (int a = p_hc.hc_art_elig_age(t); a < p_op.p_idx_fertility_first; ++a) {
        for (int hd = 0; hd < ss_c.hc1DS; ++hd) {
          if (hd > p_hc.hc_art_elig_cd4(a, t)) {
            if (a < ss_c.hc2_agestart) {
              state_next.children.hc_art_need_init(hd, cat, a, s) += state_next.children.hc1_hiv_pop(hd, cat, a, s);
            } else if (hd < ss_c.hc2DS) {
              state_next.children.hc_art_need_init(hd, cat, a, s) += state_next.children.hc2_hiv_pop(hd, cat, a - ss_c.hc2_agestart, s);
            }
          }
        } // end ss_c.hc1DS
      } // end a
    } // end hcTT
  } // end ss.NS
}


template<typename ModelVariant, typename real_type>
void need_for_cotrim(int t,
                     const Parameters<ModelVariant, real_type> &pars,
                     const State<ModelVariant, real_type> &state_curr,
                     State<ModelVariant, real_type> &state_next,
                     IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "need_for_cotrim can only be called for model variants where run_child_model is true");
  constexpr auto ss_d = StateSpace<ModelVariant>().dp;
  constexpr auto ss_c = StateSpace<ModelVariant>().children;
  const auto& p_hc = pars.children.children;
  const auto& p_op = pars.options;
  auto& n_hc = state_next.children;
  auto& i_hc = intermediate.children;
  const auto& c_hc = state_curr.children;

  //potentially this needs to be changed to include children started on ART
  // Births from the last 18 months are eligible
  n_hc.ctx_need = n_hc.hiv_births * 1.5;

  //All children 1.5-4 eligible
  for (int s = 0; s < ss_d.NS; ++s) {
    for (int a = 1; a < ss_c.hc2_agestart; ++a) {
      for (int cat = 0; cat < ss_c.hcTT; ++cat) {
        for (int hd = 0; hd < ss_c.hc1DS; ++hd) {
          if(a == 1){
            n_hc.ctx_need  += n_hc.hc1_hiv_pop(hd, cat, a, s) * 0.5;
          }else{
            n_hc.ctx_need  += n_hc.hc1_hiv_pop(hd, cat, a, s);
          }
        } // end ss_c.hc1DS
      } // end a
    } // end hcTT
  } // end ss_d.NS

  //All ART eligible children ages 5-14 eligible
  // all children under a certain age eligible for ART
  if(p_hc.hc_art_elig_age(t) > ss_c.hc2_agestart){
    for (int s = 0; s < ss_d.NS; ++s) {
      for (int a = 0; a < p_hc.hc_art_elig_age(t); ++a) {
        for (int cat = 0; cat < ss_c.hcTT; ++cat) {
          for (int hd = 0; hd < ss_c.hc1DS; ++hd) {
            if (a < ss_c.hc2_agestart) {
              n_hc.ctx_need  += n_hc.hc1_hiv_pop(hd, cat, a, s);
            } else if (hd < ss_c.hc2DS) {
              n_hc.ctx_need  += n_hc.hc2_hiv_pop(hd, cat, a - ss_c.hc2_agestart, s);
            }
          } // end ss_c.hc1DS
        } // end a
      } // end hcTT
    } // end ss_d.NS
  }

  // all children under a certain CD4 eligible for ART
  //Spectrum uses a lagged population and eligibility for children over five which I think it a mistake
  for (int s = 0; s < ss_d.NS; ++s) {
    for (int cat = 0; cat < ss_c.hcTT; ++cat) {
      for (int a = (ss_c.hc2_agestart); a < p_op.p_idx_fertility_first; ++a) {
        for (int hd = 0; hd < ss_c.hc2DS; ++hd) {
          if (hd > p_hc.hc_art_elig_cd4(a, t-1)) {
               n_hc.ctx_need += c_hc.hc2_hiv_pop(hd, cat, a - ss_c.hc2_agestart, s);
            }
        } // end ss_c.hc1DS
      } // end a
    } // end hcTT
  } // end ss.NS

}

template<typename ModelVariant, typename real_type>
void get_cotrim_effect(int t,
                          const Parameters<ModelVariant, real_type> &pars,
                          const State<ModelVariant, real_type> &state_curr,
                          State<ModelVariant, real_type> &state_next,
                          IntermediateData<ModelVariant, real_type> &intermediate,
                          int art_flag) {
  static_assert(ModelVariant::run_child_model,
                "get_cotrim_effect can only be called for model variants where run_child_model is true");
  const auto& p_hc = pars.children.children;
  auto& n_hc = state_next.children;

  //note this is just for off art, need to also do for on art
  if (p_hc.ctx_val_is_percent(t)) {
    n_hc.ctx_mean(art_flag) = (1 - p_hc.ctx_effect(art_flag)*  p_hc.ctx_val(t)) ;
  } else {
    if (n_hc.ctx_need > 0) {
     n_hc.ctx_mean(art_flag) = (1 - p_hc.ctx_effect(art_flag) * (p_hc.ctx_val(t) / n_hc.ctx_need));
    }
  }

}

template<typename ModelVariant, typename real_type>
void cd4_mortality(int t,
                   const Parameters<ModelVariant, real_type> &pars,
                   const State<ModelVariant, real_type> &state_curr,
                   State<ModelVariant, real_type> &state_next,
                   IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "cd4_mortality can only be called for model variants where run_child_model is true");
  constexpr auto ss_d = StateSpace<ModelVariant>().dp;
  constexpr auto ss_c = StateSpace<ModelVariant>().children;
  const auto& p_hc = pars.children.children;
  const auto& p_op = pars.options;
  auto& n_hc = state_next.children;
  auto& i_hc = intermediate.children;
  auto art_flag = 0;

  internal::get_cotrim_effect(t, pars, state_curr, state_next, intermediate, art_flag);
    for (int s = 0; s < ss_d.NS; ++s) {

    for (int a = 0; a < ss_c.hc2_agestart; ++a) {
      for (int cat = 0; cat < ss_c.hcTT; ++cat) {
        for (int hd = 0; hd < ss_c.hc1DS; ++hd) {
          auto hiv_deaths_strat = n_hc.ctx_mean(art_flag) * n_hc.hc1_hiv_pop(hd, cat, a, s) * p_hc.hc1_cd4_mort(hd, cat, a);
          i_hc.hc_posthivmort(hd, cat, a, s) = n_hc.hc1_hiv_pop(hd, cat, a, s) - hiv_deaths_strat;
        }
      }
    }
  }

  for (int s = 0; s < ss_d.NS; ++s) {
    for (int a = ss_c.hc2_agestart; a < p_op.p_idx_fertility_first; ++a) {
      for (int cat = 0; cat < ss_c.hcTT; ++cat) {
        for (int hd = 0; hd < ss_c.hc2DS; ++hd) {
          auto hiv_deaths_strat = n_hc.ctx_mean(art_flag) * n_hc.hc2_hiv_pop(hd, cat, a - ss_c.hc2_agestart, s) *
                                  p_hc.hc2_cd4_mort(hd, cat, a - ss_c.hc2_agestart);
          i_hc.hc_posthivmort(hd, cat, a, s) = n_hc.hc2_hiv_pop(hd, cat, a - ss_c.hc2_agestart, s) -
                                               hiv_deaths_strat;
        }
      }
    }
  }

  // progress through CD4 categories
  for (int s = 0; s < ss_d.NS; ++s) {
    for (int hd = 1; hd < ss_c.hc1DS; ++hd) {
      for (int a = 0; a < ss_c.hc2_agestart; ++a) {
        for (int cat = 0; cat < ss_c.hcTT; ++cat) {
          const auto& coarse_hc1_cd4_prog = p_hc.hc1_cd4_prog(hd - 1, p_hc.hc_age_coarse_cd4(a), s);
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
  for (int s = 0; s < ss_d.NS; ++s) {
    for (int hd = 1; hd < ss_c.hc2DS; ++hd) {
      for (int a = ss_c.hc2_agestart; a < p_op.p_idx_fertility_first; ++a) {
        for (int cat = 0; cat < ss_c.hcTT; ++cat) {
          auto cd4_grad = p_hc.hc2_cd4_prog(hd - 1, 0, s) *
                          (i_hc.hc_posthivmort(hd - 1, cat, a, s) + n_hc.hc2_hiv_pop(hd - 1, cat, a - ss_c.hc2_agestart, s)) /
                          2.0;
          i_hc.hc_grad(hd - 1, cat, a, s) -= cd4_grad; // moving to next cd4 category
          i_hc.hc_grad(hd, cat, a, s) += cd4_grad; // moving into this cd4 category
        }
      }
    }
  }
}

template<typename ModelVariant, typename real_type>
void run_child_hiv_mort(int t,
                        const Parameters<ModelVariant, real_type> &pars,
                        const State<ModelVariant, real_type> &state_curr,
                        State<ModelVariant, real_type> &state_next,
                        IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "run_child_hiv_mort can only be called for model variants where run_child_model is true");
  constexpr auto ss_d = StateSpace<ModelVariant>().dp;
  constexpr auto ss_c = StateSpace<ModelVariant>().children;
  const auto& p_hc = pars.children.children;
  const auto& p_op = pars.options;
  auto& n_hc = state_next.children;
  auto& i_hc = intermediate.children;
  auto art_flag = 0;

  internal::get_cotrim_effect(t, pars, state_curr, state_next, intermediate, art_flag);
  for (int s = 0; s < ss_d.NS; ++s) {
    for (int a = 0; a < ss_c.hc2_agestart; ++a) {
      for (int cat = 0; cat < ss_c.hcTT; ++cat) {
        for (int hd = 0; hd < ss_c.hc1DS; ++hd) {
          auto cd4_deaths_grad = n_hc.ctx_mean(art_flag) * n_hc.hc1_hiv_pop(hd, cat, a, s) *
                                 p_hc.hc1_cd4_mort(hd, cat, a);
          i_hc.hc_grad(hd, cat, a, s) -= cd4_deaths_grad;
          n_hc.hc1_noart_aids_deaths(hd, cat, a, s) += cd4_deaths_grad;
        }
      }
    }
  }

  for (int s = 0; s < ss_d.NS; ++s) {
    for (int a = ss_c.hc2_agestart; a < p_op.p_idx_fertility_first; ++a) {
      for (int cat = 0; cat < ss_c.hcTT; ++cat) {
        for (int hd = 0; hd < ss_c.hc2DS; ++hd) {
          auto cd4_mort_grad = n_hc.ctx_mean(art_flag) *
                               n_hc.hc2_hiv_pop(hd, cat, a - ss_c.hc2_agestart, s) *
                               p_hc.hc2_cd4_mort(hd, cat, a - ss_c.hc2_agestart);
          i_hc.hc_grad(hd, cat, a, s) -= cd4_mort_grad;
          n_hc.hc2_noart_aids_deaths(hd, cat, a - ss_c.hc2_agestart, s) += cd4_mort_grad;
        }
      }
    }
  }
}

template<typename ModelVariant, typename real_type>
void add_child_grad(int t,
                    const Parameters<ModelVariant, real_type> &pars,
                    const State<ModelVariant, real_type> &state_curr,
                    State<ModelVariant, real_type> &state_next,
                    IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "add_child_grad can only be called for model variants where run_child_model is true");
  constexpr auto ss_d = StateSpace<ModelVariant>().dp;
  constexpr auto ss_c = StateSpace<ModelVariant>().children;
  auto& n_hc = state_next.children;
  auto& i_hc = intermediate.children;

  // add on transitions
  for (int s = 0; s < ss_d.NS; ++s) {
    for (int a = 0; a < ss_c.hc2_agestart; ++a) {
      for (int cat = 0; cat < ss_c.hcTT; ++cat) {
        for (int hd = 0; hd < ss_c.hc1DS; ++hd) {
          n_hc.hc1_hiv_pop(hd, cat, a, s) += i_hc.hc_grad(hd, cat, a, s);
        } // end hc1DS
      } // end cat
    } // end a
  } // end s

  for (int s = 0; s < ss_d.NS; ++s) {
    for (int a = ss_c.hc2_agestart; a < pars.options.p_idx_fertility_first; ++a) {
      for (int cat = 0; cat < ss_c.hcTT; ++cat) {
        for (int hd = 0; hd < ss_c.hc2DS; ++hd) {
          n_hc.hc2_hiv_pop(hd, cat, a - ss_c.hc2_agestart, s) += i_hc.hc_grad(hd, cat, a, s);
        } //end hc2DS
      } //end cat
    } //end a
  } // end s
}

template<typename ModelVariant, typename real_type>
void eligible_for_treatment(int t,
                            const Parameters<ModelVariant, real_type> &pars,
                            const State<ModelVariant, real_type> &state_curr,
                            State<ModelVariant, real_type> &state_next,
                            IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "eligible_for_treatment can only be called for model variants where run_child_model is true");
  constexpr auto ss_d = StateSpace<ModelVariant>().dp;
  constexpr auto ss_c = StateSpace<ModelVariant>().children;
  const auto& p_op = pars.options;
  const auto& p_hc = pars.children.children;

  auto& n_hc = state_next.children;
  auto& i_hc = intermediate.children;

  internal::art_eligibility_by_age(t, pars, state_curr, state_next, intermediate);
  internal::art_eligibility_by_cd4(t, pars, state_curr, state_next, intermediate);

  for (int s = 0; s < ss_d.NS; ++s) {
    for (int cat = 0; cat < ss_c.hcTT; ++cat) {
      for (int a = 0; a < p_op.p_idx_fertility_first; ++a) {
        for (int hd = 0; hd < ss_c.hc1DS; ++hd) {
          i_hc.eligible(hd, a, s) += n_hc.hc_art_need_init(hd, cat, a, s);
        } // end ss_c.hc1DS
      } // end a
    } // end hcTT
  } // end ss_d.NS
}

template<typename ModelVariant, typename real_type>
void on_art_mortality(int t,
                     const Parameters<ModelVariant, real_type> &pars,
                     const State<ModelVariant, real_type> &state_curr,
                     State<ModelVariant, real_type> &state_next,
                     IntermediateData<ModelVariant, real_type> &intermediate,
                     int t_art_idx) {
  static_assert(ModelVariant::run_child_model,
                "on_art_mortality can only be called for model variants where run_child_model is true");
  constexpr auto ss_d = StateSpace<ModelVariant>().dp;
  constexpr auto ss_c = StateSpace<ModelVariant>().children;
  const auto& p_hc = pars.children.children;
  const auto& p_op = pars.options;
  auto& n_hc = state_next.children;
  auto& i_hc = intermediate.children;

  for (int s = 0; s < ss_d.NS; ++s) {
    for (int a = 0; a < p_op.p_idx_fertility_first; ++a) {
      for (int hd = 0; hd < ss_c.hc1DS; ++hd) {
        i_hc.hc_death_rate = 0.0;
        i_hc.hc_art_grad(t_art_idx, hd, a, s) = 0.0;
        bool is_hc2_art_pop = n_hc.hc2_art_pop(t_art_idx, hd, a - ss_c.hc2_agestart, s) > 0;
        if (t_art_idx == 0) {
          if (a < ss_c.hc2_agestart) {
            i_hc.hc_death_rate = p_hc.hc_art_mort_rr(t_art_idx, a, t) *
                                 (p_hc.hc1_art_mort(hd, 0, a) + p_hc.hc1_art_mort(hd, 1, a)) /
                                 2.0;
          } else if (hd < ss_c.hc2DS && is_hc2_art_pop) {
            i_hc.hc_death_rate = p_hc.hc_art_mort_rr(t_art_idx, a, t) *
                                 (p_hc.hc2_art_mort(hd, 0, a - ss_c.hc2_agestart) + p_hc.hc2_art_mort(hd, 1, a - ss_c.hc2_agestart)) /
                                 2.0;
          }
        } else {
          if (a < ss_c.hc2_agestart) {
            i_hc.hc_death_rate = p_hc.hc_art_mort_rr(t_art_idx, a, t) * p_hc.hc1_art_mort(hd, 2, a);
          } else if (hd < ss_c.hc2DS && is_hc2_art_pop) {
            i_hc.hc_death_rate = p_hc.hc_art_mort_rr(t_art_idx, a, t) * p_hc.hc2_art_mort(hd, 2, a-ss_c.hc2_agestart);
          }
        }

        bool any_hc1_art_deaths = i_hc.hc_death_rate * n_hc.hc1_art_pop(t_art_idx, hd, a, s) >= 0;
        bool any_hc2_art_deaths = i_hc.hc_death_rate * n_hc.hc2_art_pop(t_art_idx, hd, a - ss_c.hc2_agestart, s) >= 0;
        if (a < ss_c.hc2_agestart && any_hc1_art_deaths) {
          i_hc.hc_art_grad(t_art_idx,hd, a, s) -= i_hc.hc_death_rate * n_hc.hc1_art_pop(t_art_idx, hd, a, s);
          n_hc.hc1_art_pop(t_art_idx, hd,  a, s) += i_hc.hc_art_grad(t_art_idx, hd, a, s);
          n_hc.hc1_art_aids_deaths(t_art_idx,hd, a, s) -= i_hc.hc_art_grad(t_art_idx,hd, a, s);
        } else if (hd < ss_c.hc2DS && any_hc2_art_deaths) {
          i_hc.hc_art_grad(t_art_idx, hd, a, s) -= i_hc.hc_death_rate *
                                                   n_hc.hc2_art_pop(t_art_idx, hd, a - ss_c.hc2_agestart, s);
          n_hc.hc2_art_pop(t_art_idx, hd,  a-ss_c.hc2_agestart, s) += i_hc.hc_art_grad(t_art_idx, hd, a, s);
          n_hc.hc2_art_aids_deaths(t_art_idx, hd, a-ss_c.hc2_agestart, s) -= i_hc.hc_art_grad(t_art_idx, hd, a, s);
        }
      } // end a
    } // end ss_c.hc1DS
  } // end ss_d.NS
}

template<typename ModelVariant, typename real_type>
void deaths_this_year(int t,
                      const Parameters<ModelVariant, real_type> &pars,
                      const State<ModelVariant, real_type> &state_curr,
                      State<ModelVariant, real_type> &state_next,
                      IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "deaths_this_year can only be called for model variants where run_child_model is true");
  constexpr auto ss_h = StateSpace<ModelVariant>().hiv;
  constexpr auto ss_d = StateSpace<ModelVariant>().dp;
  constexpr auto ss_c = StateSpace<ModelVariant>().children;
  const auto& p_hc = pars.children.children;
  const auto& p_op = pars.options;
  auto& n_hc = state_next.children;
  auto& i_hc = intermediate.children;

  for (int dur = 0; dur < ss_h.hTS; ++dur) {
    for (int s = 0; s < ss_d.NS; ++s) {
      for (int hd = 0; hd < ss_c.hc1DS; ++hd) {
        for (int a = 0; a < p_op.p_idx_fertility_first; ++a) {
          if (a < ss_c.hc2_agestart) {
            i_hc.hc_art_deaths(p_hc.hc_age_coarse(a)) += n_hc.hc1_art_aids_deaths(dur, hd, a, s);
          } else if (hd < ss_c.hc2DS) {
            i_hc.hc_art_deaths(p_hc.hc_age_coarse(a)) += n_hc.hc2_art_aids_deaths(dur, hd, a - ss_c.hc2_agestart, s);
          }
        } // end a
      } // end ss_c.hc1DS
    } // end ss_d.NS
  } // end dur
  i_hc.hc_art_deaths(0) = i_hc.hc_art_deaths(1) + i_hc.hc_art_deaths(2) + i_hc.hc_art_deaths(3);
}

template<typename ModelVariant, typename real_type>
void progress_time_on_art(int t,
                          const Parameters<ModelVariant, real_type> &pars,
                          const State<ModelVariant, real_type> &state_curr,
                          State<ModelVariant, real_type> &state_next,
                          IntermediateData<ModelVariant, real_type> &intermediate,
                          int curr_t_idx, int end_t_idx) {
  static_assert(ModelVariant::run_child_model,
                "progress_time_on_art can only be called for model variants where run_child_model is true");
  constexpr auto ss_d = StateSpace<ModelVariant>().dp;
  constexpr auto ss_c = StateSpace<ModelVariant>().children;
  const auto& p_op = pars.options;
  auto& n_hc = state_next.children;

  // Progress ART to the correct time on ART
  for (int hd = 0; hd < ss_c.hc1DS; ++hd) {
    for (int a = 0; a < p_op.p_idx_fertility_first; ++a) {
      for (int s = 0; s < ss_d.NS; ++s) {
        if (a < ss_c.hc2_agestart) {
          if (n_hc.hc1_art_pop(curr_t_idx, hd, a, s) > 0) {
            n_hc.hc1_art_pop(end_t_idx, hd, a, s) += n_hc.hc1_art_pop(curr_t_idx, hd, a, s);
            n_hc.hc1_art_pop(curr_t_idx, hd, a, s) -= n_hc.hc1_art_pop(curr_t_idx, hd, a, s);
          }
        } else if (hd < ss_c.hc2DS) {
          n_hc.hc2_art_pop(end_t_idx, hd, a - ss_c.hc2_agestart, s) += n_hc.hc2_art_pop(curr_t_idx, hd, a - ss_c.hc2_agestart, s);
          n_hc.hc2_art_pop(curr_t_idx, hd, a - ss_c.hc2_agestart, s) -= n_hc.hc2_art_pop(curr_t_idx, hd, a - ss_c.hc2_agestart, s);
        }
      } // end ss_d.NS
    } // end a
  } // end ss_c.hc1DS
}

template<typename ModelVariant, typename real_type>
void calc_total_and_unmet_art_need(int t,
                                   const Parameters<ModelVariant, real_type> &pars,
                                   const State<ModelVariant, real_type> &state_curr,
                                   State<ModelVariant, real_type> &state_next,
                                   IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "calc_total_and_unmet_art_need can only be called for model variants where run_child_model is true");
  constexpr auto ss_d = StateSpace<ModelVariant>().dp;
  constexpr auto ss_h = StateSpace<ModelVariant>().hiv;
  constexpr auto ss_c = StateSpace<ModelVariant>().children;
  const auto& p_hc = pars.children.children;
  const auto& p_op = pars.options;
  auto& n_hc = state_next.children;
  auto& i_hc = intermediate.children;

  internal::eligible_for_treatment(t, pars, state_curr, state_next, intermediate);

  for (int s = 0; s < ss_d.NS; ++s) {
    for (int a = 0; a < p_op.p_idx_fertility_first; ++a) {
      for (int hd = 0; hd < ss_c.hc1DS; ++hd) {
        i_hc.unmet_need(p_hc.hc_age_coarse(a)) += i_hc.eligible(hd, a, s) ;
      } // end ss_c.hc1DS
    } // end a
  } // end ss_d.NS

  for (int s = 0; s < ss_d.NS; ++s) {
    for (int a = 0; a < p_op.p_idx_fertility_first; ++a) {
      for (int hd = 0; hd < ss_c.hc1DS; ++hd) {
        for (int dur = 0; dur < ss_h.hTS; ++dur) {
          if (a < ss_c.hc2_agestart) {
            i_hc.on_art(p_hc.hc_age_coarse(a)) += n_hc.hc1_art_pop(dur, hd, a, s);
          } else if (hd < (ss_c.hc2DS)) {
            i_hc.on_art(p_hc.hc_age_coarse(a)) += n_hc.hc2_art_pop(dur, hd, a - ss_c.hc2_agestart, s);
          }
        }
      }
    }
  }

  for (int ag = 1; ag < 4; ++ag) {
    i_hc.on_art(0) += i_hc.on_art(ag);
    i_hc.unmet_need(0) += i_hc.unmet_need(ag);
    i_hc.total_need(ag) += i_hc.on_art(ag) + i_hc.unmet_need(ag) + i_hc.hc_art_deaths(ag);
  } // end ag
  i_hc.total_need(0) = i_hc.on_art(0) + i_hc.unmet_need(0) + i_hc.hc_art_deaths(0);
}

template<typename ModelVariant, typename real_type>
void age_specific_art_last_year(int t,
                                 const Parameters<ModelVariant, real_type> &pars,
                                 const State<ModelVariant, real_type> &state_curr,
                                 State<ModelVariant, real_type> &state_next,
                                 IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "age_specific_art_last_year can only be called for model variants where run_child_model is true");
  constexpr auto ss_d = StateSpace<ModelVariant>().dp;
  constexpr auto ss_h = StateSpace<ModelVariant>().hiv;
  constexpr auto ss_c = StateSpace<ModelVariant>().children;
  const auto& p_hc = pars.children.children;
  const auto& p_op = pars.options;
  const auto& c_hc = state_curr.children;
  auto& i_hc = intermediate.children;

  if (p_hc.hc_art_is_age_spec(t - 1)) {
    for (int ag = 1; ag < 4; ++ag) {
      i_hc.total_art_last_year(ag) = p_hc.hc_art_val(ag, t - 1);
    } // end ag
  } else {
    for (int s = 0; s < ss_d.NS; ++s) {
      for (int a = 0; a < p_op.p_idx_fertility_first; ++a) {
        for (int hd = 0; hd < ss_c.hc1DS; ++hd) {
          for (int dur = 0; dur < ss_h.hTS; ++dur) {
            if (a < ss_c.hc2_agestart) {
              i_hc.total_art_last_year(p_hc.hc_age_coarse(a)) += c_hc.hc1_art_pop(dur, hd, a, s);
            } else if (hd < ss_c.hc2DS) {
              i_hc.total_art_last_year(p_hc.hc_age_coarse(a)) += c_hc.hc2_art_pop(dur, hd, a - ss_c.hc2_agestart, s);
            }
          }
        }
      }
    }
    i_hc.total_art_last_year(0) = i_hc.total_art_last_year(1) +
                                  i_hc.total_art_last_year(2) +
                                  i_hc.total_art_last_year(3);

    for (int ag = 1; ag < 4; ++ag) {
      i_hc.total_art_last_year(ag) = p_hc.hc_art_val(0, t - 1) *
                                     i_hc.total_art_last_year(ag) / i_hc.total_art_last_year(0);
      if (p_hc.hc_art_isperc(t - 1)) {
        i_hc.total_art_last_year(ag) *= i_hc.total_need(0) + i_hc.hc_art_deaths(ag);
      }
    } // end ag
  }
}

template<typename ModelVariant, typename real_type>
void art_last_year(int t,
                   const Parameters<ModelVariant, real_type> &pars,
                   const State<ModelVariant, real_type> &state_curr,
                   State<ModelVariant, real_type> &state_next,
                   IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "art_last_year can only be called for model variants where run_child_model is true");
  const auto& p_hc = pars.children.children;
  auto& i_hc = intermediate.children;

  if (p_hc.hc_art_is_age_spec(t)) {
    // If the present time step is age specific, we need to calculate what last years age spec
    // breakdown would have been

    // Age specific ART will always be entered as a number
    internal::age_specific_art_last_year(t, pars, state_curr, state_next, intermediate);
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
}

template<typename ModelVariant, typename real_type>
void art_this_year(int t,
                   const Parameters<ModelVariant, real_type> &pars,
                   const State<ModelVariant, real_type> &state_curr,
                   State<ModelVariant, real_type> &state_next,
                   IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "art_this_year can only be called for model variants where run_child_model is true");
  const auto& p_hc = pars.children.children;
  auto& i_hc = intermediate.children;

  for (int ag = 0; ag < 4; ++ag) {
    i_hc.total_art_this_year(ag) = p_hc.hc_art_val(ag, t);
    if (p_hc.hc_art_isperc(t)) {
      i_hc.total_art_this_year(ag) *= i_hc.total_need(ag);
    }
  } // end ag
}

template<typename ModelVariant, typename real_type>
void calc_art_initiates(int t,
                        const Parameters<ModelVariant, real_type> &pars,
                        const State<ModelVariant, real_type> &state_curr,
                        State<ModelVariant, real_type> &state_next,
                        IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "calc_art_initiates can only be called for model variants where run_child_model is true");
  const auto& p_hc = pars.children.children;
  auto& n_hc = state_next.children;
  auto& i_hc = intermediate.children;

  internal::on_art_mortality(t, pars, state_curr, state_next, intermediate, 0);
  // progress art initates from 0-6 months on art to 6 to 12 mo
  internal::progress_time_on_art(t, pars, state_curr, state_next, intermediate, 0, 1);
  internal::on_art_mortality(t, pars, state_curr, state_next, intermediate, 2);
  internal::deaths_this_year(t, pars, state_curr, state_next, intermediate);
  internal::calc_total_and_unmet_art_need(t, pars, state_curr, state_next, intermediate);
  internal::art_last_year(t, pars, state_curr, state_next, intermediate);
  internal::art_this_year(t, pars, state_curr, state_next, intermediate);

  i_hc.retained = 1 - p_hc.hc_art_ltfu(t);
  for (int ag = 0; ag < 4; ++ag) {
    auto average_art_by_year = (i_hc.total_art_last_year(ag) + i_hc.total_art_this_year(ag)) /
                               2.0;
    n_hc.hc_art_init(ag) = std::max(average_art_by_year - i_hc.on_art(ag) * i_hc.retained, 0.0);
    n_hc.hc_art_init(ag) = std::min(n_hc.hc_art_init(ag),
                                    i_hc.unmet_need(ag) + i_hc.on_art(ag) * p_hc.hc_art_ltfu(t));
  } // end ag

}

template<typename ModelVariant, typename real_type>
void art_ltfu(int t,
              const Parameters<ModelVariant, real_type> &pars,
              const State<ModelVariant, real_type> &state_curr,
              State<ModelVariant, real_type> &state_next,
              IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "art_ltfu can only be called for model variants where run_child_model is true");
  constexpr auto ss_d = StateSpace<ModelVariant>().dp;
  constexpr auto ss_c = StateSpace<ModelVariant>().children;
  const auto& p_hc = pars.children.children;
  const auto& p_op = pars.options;
  auto& n_hc = state_next.children;
  auto& i_hc = intermediate.children;

  for (int s = 0; s < ss_d.NS; ++s) {
    for (int a = 0; a < p_op.p_idx_fertility_first; ++a) {
      for (int hd = 0; hd < ss_c.hc1DS; ++hd) {
        for (int cat = 0; cat < ss_c.hcTT; ++cat) {
          if (a < ss_c.hc2_agestart) {
            i_hc.hc_hiv_total(hd, a, s) += n_hc.hc1_hiv_pop(hd, cat, a, s);
          } else if (hd < ss_c.hc2DS) {
            i_hc.hc_hiv_total(hd, a, s) += n_hc.hc2_hiv_pop(hd, cat, a - ss_c.hc2_agestart, s);
          }
        } // end ss_c.hcTT
      } // end ss_c.hc1DS
    } // end a
  } // end ss_d.NS

  for (int s = 0; s <ss_d.NS; ++s) {
    for (int a = 0; a < p_op.p_idx_fertility_first; ++a) {
      for (int hd = 0; hd < ss_c.hc1DS; ++hd) {
        for (int cat = 0; cat < ss_c.hcTT; ++cat) {
          if (a < ss_c.hc2_agestart) {
            i_hc.hc_hiv_dist(hd, cat, a, s) += n_hc.hc1_hiv_pop(hd, cat, a, s) / i_hc.hc_hiv_total(hd, a, s);
          } else if (hd < ss_c.hc2DS) {
            i_hc.hc_hiv_dist(hd, cat, a, s) += n_hc.hc2_hiv_pop(hd, cat, a - ss_c.hc2_agestart, s) / i_hc.hc_hiv_total(hd, a, s);
          }
        } // end ss_c.hcTT
      } // end ss_c.hc1DS
    } // end a
  } // end ss_d.NS

  for (int s = 0; s < ss_d.NS; ++s) {
    for (int a = 0; a < p_op.p_idx_fertility_first; ++a) {
      for (int hd = 0; hd < ss_c.hc1DS; ++hd) {
        for (int cat = 0; cat < ss_c.hcTT; ++cat) {
          if (a < ss_c.hc2_agestart) {
            auto ltfu_grad = (n_hc.hc1_art_pop(2, hd, a, s) + n_hc.hc1_art_pop(0, hd, a, s)) *
                             p_hc.hc_art_ltfu(t);
            if (i_hc.hc_hiv_total(hd, a, s) > 0) {
              i_hc.art_ltfu_grad(hd, cat, a, s) += ltfu_grad * i_hc.hc_hiv_dist(hd, cat, a, s);
            } else {
              i_hc.art_ltfu_grad(hd, cat, a, s) += ltfu_grad * 0.25;
            }
          } else if (hd < ss_c.hc2DS) {
            auto ltfu_grad = (n_hc.hc2_art_pop(2, hd, a - ss_c.hc2_agestart, s) + n_hc.hc2_art_pop(0, hd, a - ss_c.hc2_agestart, s)) *
                             p_hc.hc_art_ltfu(t);
            if (i_hc.hc_hiv_total(hd, a, s) > 0) {
              i_hc.art_ltfu_grad(hd, cat, a, s) += ltfu_grad * i_hc.hc_hiv_dist(hd, cat, a, s);
            } else {
              i_hc.art_ltfu_grad(hd, cat, a, s) += ltfu_grad * 0.25;
            }
          }
        } // end ss_c.hcTT
      } // end ss_c.hc1DS
    } // end a
  } // end ss_d.NS
}

template<typename ModelVariant, typename real_type>
void apply_ltfu_to_hivpop(int t,
                          const Parameters<ModelVariant, real_type> &pars,
                          const State<ModelVariant, real_type> &state_curr,
                          State<ModelVariant, real_type> &state_next,
                          IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "apply_ltfu_to_hivpop can only be called for model variants where run_child_model is true");
  constexpr auto ss_d = StateSpace<ModelVariant>().dp;
  constexpr auto ss_c = StateSpace<ModelVariant>().children;
  const auto& p_op = pars.options;
  auto& n_hc = state_next.children;
  auto& i_hc = intermediate.children;

  for (int s = 0; s < ss_d.NS; ++s) {
    for (int a = 0; a < p_op.p_idx_fertility_first; ++a) {
      for (int hd = 0; hd < ss_c.hc1DS; ++hd) {
        for (int cat = 0; cat < ss_c.hcTT; ++cat) {
          if (a < ss_c.hc2_agestart) {
            n_hc.hc1_hiv_pop(hd, cat, a, s) += i_hc.art_ltfu_grad(hd, cat, a, s);
          } else if (hd < ss_c.hc2DS) {
            n_hc.hc2_hiv_pop(hd, cat, a - ss_c.hc2_agestart, s) += i_hc.art_ltfu_grad(hd, cat, a, s);
          }
        } // end ss_c.hcTT
      } // end ss_c.hc1DS
    } // end a
  } // end ss_d.NS
}

template<typename ModelVariant, typename real_type>
void apply_ltfu_to_artpop(int t,
                          const Parameters<ModelVariant, real_type> &pars,
                          const State<ModelVariant, real_type> &state_curr,
                          State<ModelVariant, real_type> &state_next,
                          IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "apply_ltfu_to_artpop can only be called for model variants where run_child_model is true");
  constexpr auto ss_d = StateSpace<ModelVariant>().dp;
  constexpr auto ss_c = StateSpace<ModelVariant>().children;
  const auto& p_op = pars.options;
  auto& n_hc = state_next.children;
  auto& i_hc = intermediate.children;

  for (int s = 0; s < ss_d.NS; ++s) {
    for (int a = 0; a < p_op.p_idx_fertility_first; ++a) {
      for (int hd = 0; hd < ss_c.hc1DS; ++hd) {
        for (int cat = 0; cat < ss_c.hcTT; ++cat) {
          if (a < ss_c.hc2_agestart) {
            // hard coded two as this will only occur among children that are on ART more than a year
            n_hc.hc1_art_pop(2, hd, a, s) -= i_hc.art_ltfu_grad(hd, cat, a, s);
          } else if (hd < ss_c.hc2DS) {
            n_hc.hc2_art_pop(2, hd, a - ss_c.hc2_agestart, s) -= i_hc.art_ltfu_grad(hd, cat, a, s);
          }
        } // end ss_c.hcTT
      } // end ss_c.hc1DS
    } // end a
  } // end ss_d.NS
}

template<typename ModelVariant, typename real_type>
void art_initiation_by_age(int t,
                           const Parameters<ModelVariant, real_type> &pars,
                           const State<ModelVariant, real_type> &state_curr,
                           State<ModelVariant, real_type> &state_next,
                           IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "art_initiation_by_age can only be called for model variants where run_child_model is true");
  constexpr auto ss_d = StateSpace<ModelVariant>().dp;
  constexpr auto ss_c = StateSpace<ModelVariant>().children;
  const auto& p_hc = pars.children.children;
  const auto& p_op = pars.options;
  auto& n_hc = state_next.children;
  auto& i_hc = intermediate.children;

  internal::calc_art_initiates(t, pars, state_curr, state_next, intermediate);

  if (p_hc.hc_art_is_age_spec(t)) {
    for (int s = 0; s < ss_d.NS; ++s) {
      for (int a = 0; a < p_op.p_idx_fertility_first; ++a) {
        for (int hd = 0; hd < ss_c.hc1DS; ++hd) {
          for (int cat = 0; cat < ss_c.hcTT; ++cat) {
            i_hc.hc_initByAge(p_hc.hc_age_coarse(a)) += n_hc.hc_art_need_init(hd, cat, a, s) *
                                                        p_hc.hc_art_init_dist(a, t);
          } // end hcTT
        } // end ss_c.hc1DS
      } // end a
    } // end ss_d.NS

    for (int ag = 1; ag < 4; ++ag) {
      if (i_hc.hc_initByAge(ag) == 0.0) {
        i_hc.hc_adj(ag) = 1.0;
      } else {
        i_hc.hc_adj(ag) = n_hc.hc_art_init(ag) / i_hc.hc_initByAge(ag);
      }
    }

    for (int s = 0; s <ss_d.NS; ++s) {
      for (int cat = 0; cat < ss_c.hcTT; ++cat) {
        for (int a = 0; a < p_op.p_idx_fertility_first; ++a) {
          for (int hd = 0; hd < ss_c.hc1DS; ++hd) {
            auto& coarse_hc_adj = i_hc.hc_adj(p_hc.hc_age_coarse(a));
            auto& coarse_hc_art_scalar = i_hc.hc_art_scalar(p_hc.hc_age_coarse(a));
            auto hc_art_val_sum = p_hc.hc_art_val(0, t) + p_hc.hc_art_val(1, t) +
                                  p_hc.hc_art_val(2, t) + p_hc.hc_art_val(3, t);
            auto hc_art_val_sum_last = p_hc.hc_art_val(0, t-1) + p_hc.hc_art_val(1, t-1) +
              p_hc.hc_art_val(2, t-1) + p_hc.hc_art_val(3, t-1);
            if ((hc_art_val_sum + hc_art_val_sum_last) <= 0.0) {
              coarse_hc_art_scalar = 0.0;
            } else {
              coarse_hc_art_scalar = std::min(coarse_hc_adj * p_hc.hc_art_init_dist(a, t), 1.0);
            }

            auto art_initiates = coarse_hc_art_scalar * n_hc.hc_art_need_init(hd, cat, a, s);

            if (a < ss_c.hc2_agestart) {
              n_hc.hc1_art_pop(0, hd, a, s) += art_initiates;
              n_hc.hc1_hiv_pop(hd, cat, a, s) -= art_initiates;
            } else if (hd < (ss_c.hc2DS)) {
              n_hc.hc2_art_pop(0, hd, a - ss_c.hc2_agestart, s) += art_initiates;
              n_hc.hc2_hiv_pop(hd, cat, a - ss_c.hc2_agestart, s) -=  art_initiates;
            }
          } // end ss_c.hc1DS
        } // end a
      } // end hcTT
    } // end ss_d.NS
  } else {
    for (int s = 0; s <ss_d.NS; ++s) {
      for (int a = 0; a < p_op.p_idx_fertility_first; ++a) {
        for (int hd = 0; hd < ss_c.hc1DS; ++hd) {
          for (int cat = 0; cat < ss_c.hcTT; ++cat) {
            i_hc.hc_initByAge(0) += n_hc.hc_art_need_init(hd, cat, a, s) * p_hc.hc_art_init_dist(a, t);
          } // end hcTT
        } // end ss_c.hc1DS
      } // end a
    } // end ss_d.NS

    if (i_hc.hc_initByAge(0) == 0.0) {
      i_hc.hc_adj(0) = 1.0 ;
    } else {
      i_hc.hc_adj(0) =  n_hc.hc_art_init(0) / i_hc.hc_initByAge(0);
    }

    for (int s = 0; s <ss_d.NS; ++s) {
      for (int cat = 0; cat < ss_c.hcTT; ++cat) {
        for (int a = 0; a < p_op.p_idx_fertility_first; ++a) {
          for (int hd = 0; hd < ss_c.hc1DS; ++hd) {
            auto hc_art_val_sum = p_hc.hc_art_val(0, t) + p_hc.hc_art_val(0, t - 1);
            if (hc_art_val_sum <= 0) {
              i_hc.hc_art_scalar(0) = 0.0;
            } else {
              // issue is that in 2030 this is being activated when it shouldn't be for age one
              // TLDR is that too many people are initiating who shouldn't be
              i_hc.hc_art_scalar(0) = std::min(i_hc.hc_adj(0) * p_hc.hc_art_init_dist(a, t), 1.0);
            }

            auto art_initiates = i_hc.hc_art_scalar(0) * n_hc.hc_art_need_init(hd, cat, a, s);
            if (a < ss_c.hc2_agestart) {
              n_hc.hc1_art_pop(0, hd, a, s) += art_initiates;
              n_hc.hc1_hiv_pop(hd, cat, a, s) -= art_initiates;
            } else if (hd < (ss_c.hc2DS)) {
              n_hc.hc2_art_pop(0, hd, a - ss_c.hc2_agestart, s) += art_initiates;
              n_hc.hc2_hiv_pop(hd, cat, a - ss_c.hc2_agestart, s) -= art_initiates;
            }
          } // end ss_c.hc1DS
        } // end a
      } // end hcTT
    } // end ss_d.NS
  } // end if


}

template<typename ModelVariant, typename real_type>
void fill_total_pop_outputs(int t,
                            const Parameters<ModelVariant, real_type> &pars,
                            const State<ModelVariant, real_type> &state_curr,
                            State<ModelVariant, real_type> &state_next,
                            IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "fill_total_pop_outputs can only be called for model variants where run_child_model is true");
  constexpr auto ss_d = StateSpace<ModelVariant>().dp;
  constexpr auto ss_h = StateSpace<ModelVariant>().hiv;
  constexpr auto ss_c = StateSpace<ModelVariant>().children;
  const auto& p_op = pars.options;
  auto& n_ha = state_next.hiv;
  auto& n_hc = state_next.children;

  for (int hd = 0; hd < ss_h.hDS; ++hd) {
    for (int a = 0; a < p_op.p_idx_fertility_first; ++a) {
      for (int s = 0; s < ss_d.NS; ++s) {
        for (int cat = 0; cat < ss_c.hcTT; ++cat) {
          if (a < ss_c.hc2_agestart) {
            n_ha.p_hiv_deaths(a,s) += n_hc.hc1_noart_aids_deaths(hd, cat, a, s);
          } else if (hd < ss_c.hc2DS) {
            n_ha.p_hiv_deaths(a,s) += n_hc.hc2_noart_aids_deaths(hd, cat, a - ss_c.hc2_agestart, s);
          }
        } // end hcTT

        for (int dur = 0; dur < ss_h.hTS; ++dur) {
          if (a < ss_c.hc2_agestart) {
            n_ha.p_hiv_deaths(a,s) += n_hc.hc1_art_aids_deaths(dur, hd, a, s);
          } else if (hd < ss_c.hc2DS) {
            n_ha.p_hiv_deaths(a,s) += n_hc.hc2_art_aids_deaths(dur, hd, a - ss_c.hc2_agestart, s);
          }
        } // end dur
      } // end ss_d.NS
    } // end a
  } // end hDS


  for (int a = 0; a < p_op.p_idx_fertility_first; ++a) {
  for (int s = 0; s < ss_d.NS; ++s) {
    // if(s == 1 & a > 0){
        n_ha.p_hiv_pop(a, s) += n_ha.p_infections(a, s);
       // }
      n_ha.p_hiv_pop(a, s) -= n_ha.p_hiv_deaths(a, s);
    }
  }


}

} // namespace internal

template<typename ModelVariant, typename real_type>
void run_child_model_simulation(int t,
                                const Parameters<ModelVariant, real_type> &pars,
                                const State<ModelVariant, real_type> &state_curr,
                                State<ModelVariant, real_type> &state_next,
                                internal::IntermediateData<ModelVariant, real_type> &intermediate) {
  static_assert(ModelVariant::run_child_model,
                "run_child_model_simulation can only be called for model variants where run_child_model is true");
  const auto& p_hc = pars.children.children;
  const auto& p_op = pars.options;
  auto& n_hc = state_next.children;

  internal::run_child_ageing(t, pars, state_curr, state_next, intermediate);

  if (p_hc.mat_prev_input(t)) {
    internal::run_wlhiv_births_input_mat_prev(t, pars, state_curr, state_next, intermediate);
  } else {
    internal::run_wlhiv_births(t, pars, state_curr, state_next, intermediate);
  }
  internal::adjust_hiv_births(t, pars, state_curr, state_next, intermediate);
  internal::add_infections(t, pars, state_curr, state_next, intermediate);
  internal::need_for_cotrim(t, pars, state_curr, state_next, intermediate);
  internal::cd4_mortality(t, pars, state_curr, state_next, intermediate);
  internal::run_child_hiv_mort(t, pars, state_curr, state_next, intermediate);
  internal::add_child_grad(t, pars, state_curr, state_next, intermediate);

    // assume paed art doesn't start before adult
    if (t >= p_hc.hc_art_start) {
      internal::art_ltfu(t, pars, state_curr, state_next, intermediate);
      internal::art_initiation_by_age(t, pars, state_curr, state_next, intermediate);
      // mortality among those on ART less than one year
      internal::on_art_mortality(t, pars, state_curr, state_next, intermediate, 0);
      internal::progress_time_on_art(t, pars, state_curr, state_next, intermediate, 1, 2);
      // progress 6 to 12 mo to 12 plus months#
      internal::apply_ltfu_to_hivpop(t, pars, state_curr, state_next, intermediate);
      internal::apply_ltfu_to_artpop(t, pars, state_curr, state_next, intermediate);
    }

    internal::nosocomial_infections(t, pars, state_curr, state_next, intermediate);
    internal::fill_total_pop_outputs(t, pars, state_curr, state_next, intermediate);


}

} // namespace leapfrog
