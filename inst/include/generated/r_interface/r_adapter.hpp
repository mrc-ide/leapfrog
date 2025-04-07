#pragma once

#include <Rcpp.h>

#include "../config.hpp"

namespace leapfrog {

template <typename T>
T* r_data(SEXP x) {
  static_assert(sizeof(T) == 0, "Only specializations of r_data can be used");
}

template <>
double* r_data(SEXP x) {
  return REAL(x);
}

template <>
int* r_data(SEXP x) {
  return INTEGER(x);
}

template<typename T, typename... Args>
auto parse_data(const Rcpp::List data, const std::string& key, Args... dims) {
  constexpr std::size_t rank = sizeof...(dims);
  Eigen::array<int, rank> dimensions{ static_cast<int>(dims)... };

  int length = std::accumulate(dimensions.begin(), dimensions.end(), 1, std::multiplies<int>());
  SEXP array_data = data[key];
  // In cases where the input data has project years we might not use all of it model fit
  // So we can take create a Map over a smaller slice of the data
  // As long as this is true we can be confident we're not referencing invalid memory
  if (LENGTH(array_data) < length) {
    Rcpp::stop("Invalid size of data for '%s', expected %d got %d",
               key,
               length,
               LENGTH(array_data));
  }

  return Eigen::TensorMap<Eigen::Tensor<T, rank>>(r_data<T>(array_data), static_cast<int>(dims)...);
}

template<typename real_type, MV ModelVariant>
struct DpAdapterR {
  using SS = SSMixed<ModelVariant>;
  using Config = DpConfig<real_type, ModelVariant>;

  static Config::Pars get_pars(
    const Rcpp::List &data,
    const Opts<real_type> &opts,
    const int proj_years
  ) {
    return {
      .base_pop = parse_data<real_type>(data, "basepop", SS::pAG, SS::NS),
      .survival_probability = parse_data<real_type>(data, "Sx", SS::pAG + 1, SS::NS, proj_years),
      .net_migration = parse_data<real_type>(data, "netmigr_adj", SS::pAG, SS::NS, proj_years),
      .age_specific_fertility_rate = parse_data<real_type>(data, "asfr", opts.p_fertility_age_groups, proj_years),
      .births_sex_prop = parse_data<real_type>(data, "births_sex_prop", SS::NS, proj_years)
    };
  };

  static constexpr int output_count = 3;

  static int build_output(
    Rcpp::List& ret,
    Rcpp::CharacterVector& names,
    int index,
    const Config::OutputState& state,
    const size_t& output_years
  ) {
    Rcpp::NumericVector r_p_total_pop(SS::pAG * SS::NS * output_years);
    r_p_total_pop.attr("dim") = Rcpp::IntegerVector::create(SS::pAG, SS::NS, output_years);
    std::copy_n(state.p_total_pop.data(), state.p_total_pop.size(), REAL(r_p_total_pop));
    names[index + 0] = "p_total_pop";
    ret[index + 0] = r_p_total_pop;
    Rcpp::NumericVector r_p_total_pop_natural_deaths(SS::pAG * SS::NS * output_years);
    r_p_total_pop_natural_deaths.attr("dim") = Rcpp::IntegerVector::create(SS::pAG, SS::NS, output_years);
    std::copy_n(state.p_total_pop_natural_deaths.data(), state.p_total_pop_natural_deaths.size(), REAL(r_p_total_pop_natural_deaths));
    names[index + 1] = "p_total_pop_natural_deaths";
    ret[index + 1] = r_p_total_pop_natural_deaths;
    Rcpp::NumericVector r_births(output_years);
    r_births.attr("dim") = Rcpp::IntegerVector::create(output_years);
    std::copy_n(state.births.data(), state.births.size(), REAL(r_births));
    names[index + 2] = "births";
    ret[index + 2] = r_births;
    return index + output_count;
  };
};

template<typename real_type, MV ModelVariant>
struct HaAdapterR {
  using SS = SSMixed<ModelVariant>;
  using Config = HaConfig<real_type, ModelVariant>;

  static Config::Pars get_pars(
    const Rcpp::List &data,
    const Opts<real_type> &opts,
    const int proj_years
  ) {
    return {
      .total_rate = parse_data<real_type>(data, "incidinput", proj_years),
      .relative_risk_age = parse_data<real_type>(data, "incrr_age", SS::pAG - opts.p_idx_hiv_first_adult, SS::NS, proj_years),
      .relative_risk_sex = parse_data<real_type>(data, "incrr_sex", proj_years),
      .cd4_mortality = parse_data<real_type>(data, "cd4_mort", SS::hDS, SS::hAG, SS::NS),
      .cd4_progression = parse_data<real_type>(data, "cd4_prog", SS::hDS - 1, SS::hAG, SS::NS),
      .cd4_initial_distribution = parse_data<real_type>(data, "cd4_initdist", SS::hDS, SS::hAG, SS::NS),
      .scale_cd4_mortality = Rcpp::as<int>(data["scale_cd4_mort"]),
      .idx_hm_elig = parse_data<int>(data, "artcd4elig_idx", proj_years),
      .mortality = parse_data<real_type>(data, "art_mort", SS::hTS, SS::hDS, SS::hAG, SS::NS),
      .mortality_time_rate_ratio = parse_data<real_type>(data, "artmx_timerr", SS::hTS, proj_years),
      .dropout_recover_cd4 = Rcpp::as<int>(data["art_dropout_recover_cd4"]),
      .dropout_rate = parse_data<real_type>(data, "art_dropout_rate", proj_years),
      .adults_on_art = parse_data<real_type>(data, "art15plus_num", SS::NS, proj_years),
      .adults_on_art_is_percent = parse_data<int>(data, "art15plus_isperc", SS::NS, proj_years),
      .h_art_stage_dur = parse_data<real_type>(data, "h_art_stage_dur", SS::hTS - 1),
      .initiation_mortality_weight = Rcpp::as<real_type>(data["art_alloc_mxweight"])
    };
  };

  static constexpr int output_count = 9;

  static int build_output(
    Rcpp::List& ret,
    Rcpp::CharacterVector& names,
    int index,
    const Config::OutputState& state,
    const size_t& output_years
  ) {
    Rcpp::NumericVector r_p_hiv_pop(SS::pAG * SS::NS * output_years);
    r_p_hiv_pop.attr("dim") = Rcpp::IntegerVector::create(SS::pAG, SS::NS, output_years);
    std::copy_n(state.p_hiv_pop.data(), state.p_hiv_pop.size(), REAL(r_p_hiv_pop));
    names[index + 0] = "p_hiv_pop";
    ret[index + 0] = r_p_hiv_pop;
    Rcpp::NumericVector r_p_hiv_pop_natural_deaths(SS::pAG * SS::NS * output_years);
    r_p_hiv_pop_natural_deaths.attr("dim") = Rcpp::IntegerVector::create(SS::pAG, SS::NS, output_years);
    std::copy_n(state.p_hiv_pop_natural_deaths.data(), state.p_hiv_pop_natural_deaths.size(), REAL(r_p_hiv_pop_natural_deaths));
    names[index + 1] = "p_hiv_pop_natural_deaths";
    ret[index + 1] = r_p_hiv_pop_natural_deaths;
    Rcpp::NumericVector r_h_hiv_adult(SS::hDS * SS::hAG * SS::NS * output_years);
    r_h_hiv_adult.attr("dim") = Rcpp::IntegerVector::create(SS::hDS, SS::hAG, SS::NS, output_years);
    std::copy_n(state.h_hiv_adult.data(), state.h_hiv_adult.size(), REAL(r_h_hiv_adult));
    names[index + 2] = "h_hiv_adult";
    ret[index + 2] = r_h_hiv_adult;
    Rcpp::NumericVector r_h_art_adult(SS::hTS * SS::hDS * SS::hAG * SS::NS * output_years);
    r_h_art_adult.attr("dim") = Rcpp::IntegerVector::create(SS::hTS, SS::hDS, SS::hAG, SS::NS, output_years);
    std::copy_n(state.h_art_adult.data(), state.h_art_adult.size(), REAL(r_h_art_adult));
    names[index + 3] = "h_art_adult";
    ret[index + 3] = r_h_art_adult;
    Rcpp::NumericVector r_h_hiv_deaths_no_art(SS::hDS * SS::hAG * SS::NS * output_years);
    r_h_hiv_deaths_no_art.attr("dim") = Rcpp::IntegerVector::create(SS::hDS, SS::hAG, SS::NS, output_years);
    std::copy_n(state.h_hiv_deaths_no_art.data(), state.h_hiv_deaths_no_art.size(), REAL(r_h_hiv_deaths_no_art));
    names[index + 4] = "h_hiv_deaths_no_art";
    ret[index + 4] = r_h_hiv_deaths_no_art;
    Rcpp::NumericVector r_p_infections(SS::pAG * SS::NS * output_years);
    r_p_infections.attr("dim") = Rcpp::IntegerVector::create(SS::pAG, SS::NS, output_years);
    std::copy_n(state.p_infections.data(), state.p_infections.size(), REAL(r_p_infections));
    names[index + 5] = "p_infections";
    ret[index + 5] = r_p_infections;
    Rcpp::NumericVector r_h_hiv_deaths_art(SS::hTS * SS::hDS * SS::hAG * SS::NS * output_years);
    r_h_hiv_deaths_art.attr("dim") = Rcpp::IntegerVector::create(SS::hTS, SS::hDS, SS::hAG, SS::NS, output_years);
    std::copy_n(state.h_hiv_deaths_art.data(), state.h_hiv_deaths_art.size(), REAL(r_h_hiv_deaths_art));
    names[index + 6] = "h_hiv_deaths_art";
    ret[index + 6] = r_h_hiv_deaths_art;
    Rcpp::NumericVector r_h_art_initiation(SS::hDS * SS::hAG * SS::NS * output_years);
    r_h_art_initiation.attr("dim") = Rcpp::IntegerVector::create(SS::hDS, SS::hAG, SS::NS, output_years);
    std::copy_n(state.h_art_initiation.data(), state.h_art_initiation.size(), REAL(r_h_art_initiation));
    names[index + 7] = "h_art_initiation";
    ret[index + 7] = r_h_art_initiation;
    Rcpp::NumericVector r_p_hiv_deaths(SS::pAG * SS::NS * output_years);
    r_p_hiv_deaths.attr("dim") = Rcpp::IntegerVector::create(SS::pAG, SS::NS, output_years);
    std::copy_n(state.p_hiv_deaths.data(), state.p_hiv_deaths.size(), REAL(r_p_hiv_deaths));
    names[index + 8] = "p_hiv_deaths";
    ret[index + 8] = r_p_hiv_deaths;
    return index + output_count;
  };
};


template<typename real_type, MV ModelVariant>
struct HcAdapterR {
  using SS = SSMixed<ModelVariant>;
  using Config = HcConfig<real_type, ModelVariant>;

  static Config::Pars get_pars(
    const Rcpp::List &data,
    const Opts<real_type> &opts,
    const int proj_years
  ) {
    return {
      .hc_nosocomial = parse_data<real_type>(data, "paed_incid_input", proj_years),
      .hc1_cd4_dist = parse_data<real_type>(data, "paed_cd4_dist", SS::hc2DS),
      .hc_cd4_transition = parse_data<real_type>(data, "paed_cd4_transition", SS::hc2DS, SS::hc1DS),
      .hc1_cd4_mort = parse_data<real_type>(data, "paed_cd4_mort", SS::hc1DS, SS::hcTT, SS::hc1AG),
      .hc2_cd4_mort = parse_data<real_type>(data, "adol_cd4_mort", SS::hc2DS, SS::hcTT, SS::hc2AG),
      .hc1_cd4_prog = parse_data<real_type>(data, "paed_cd4_prog", SS::hc1DS, SS::hc1AG_c, SS::NS),
      .hc2_cd4_prog = parse_data<real_type>(data, "adol_cd4_prog", SS::hc2DS, SS::hc2AG_c, SS::NS),
      .ctx_val = parse_data<real_type>(data, "ctx_val", proj_years),
      .hc_art_elig_age = parse_data<int>(data, "paed_art_elig_age", proj_years),
      .hc_art_elig_cd4 = parse_data<real_type>(data, "paed_art_elig_cd4", opts.p_idx_hiv_first_adult, proj_years),
      .hc_art_mort_rr = parse_data<real_type>(data, "mort_art_rr", SS::hTS, opts.p_idx_hiv_first_adult, proj_years),
      .hc1_art_mort = parse_data<real_type>(data, "paed_art_mort", SS::hc1DS, SS::hTS, SS::hc1AG),
      .hc2_art_mort = parse_data<real_type>(data, "adol_art_mort", SS::hc2DS, SS::hTS, SS::hc2AG),
      .hc_art_isperc = parse_data<int>(data, "artpaeds_isperc", proj_years),
      .hc_art_val = parse_data<real_type>(data, "paed_art_val", SS::hcAG_coarse, proj_years),
      .hc_art_init_dist = parse_data<real_type>(data, "init_art_dist", opts.p_idx_hiv_first_adult, proj_years),
      .adult_cd4_dist = parse_data<real_type>(data, "adult_cd4_dist", SS::hDS, SS::hc2DS),
      .fert_mult_by_age = parse_data<real_type>(data, "fert_mult_by_age", opts.p_fertility_age_groups),
      .fert_mult_off_art = parse_data<real_type>(data, "fert_mult_offart", SS::hDS),
      .fert_mult_on_art = parse_data<real_type>(data, "fert_mult_onart", opts.p_fertility_age_groups),
      .total_fertility_rate = parse_data<real_type>(data, "tfr", proj_years),
      .PMTCT = parse_data<real_type>(data, "pmtct", SS::hPS, proj_years),
      .vertical_transmission_rate = parse_data<real_type>(data, "mtct", SS::hDS + 1, SS::hVT),
      .PMTCT_transmission_rate = parse_data<real_type>(data, "pmtct_mtct", SS::hDS, SS::hPS, SS::hVT),
      .PMTCT_dropout = parse_data<real_type>(data, "pmtct_dropout", SS::hPS_dropout, proj_years),
      .PMTCT_input_is_percent = parse_data<int>(data, "pmtct_input_isperc", proj_years),
      .breastfeeding_duration_art = parse_data<real_type>(data, "bf_duration_art", SS::hBF, proj_years),
      .breastfeeding_duration_no_art = parse_data<real_type>(data, "bf_duration_no_art", SS::hBF, proj_years),
      .mat_hiv_births = parse_data<real_type>(data, "mat_hiv_births", proj_years),
      .mat_prev_input = parse_data<int>(data, "mat_prev_input", proj_years),
      .prop_lt200 = parse_data<real_type>(data, "prop_lt200", proj_years),
      .prop_gte350 = parse_data<real_type>(data, "prop_gte350", proj_years),
      .incrate = parse_data<real_type>(data, "incrate", proj_years),
      .ctx_val_is_percent = parse_data<int>(data, "ctx_val_ispercent", proj_years),
      .hc_art_is_age_spec = parse_data<int>(data, "paed_art_age_spec", proj_years),
      .hc_age_coarse = parse_data<real_type>(data, "hc_age_coarse", SS::hcAG_end),
      .abortion = parse_data<real_type>(data, "abortion", SS::hAB_ind, proj_years),
      .patients_reallocated = parse_data<real_type>(data, "patients_reallocated", proj_years),
      .hc_art_ltfu = parse_data<real_type>(data, "paed_art_ltfu", proj_years),
      .hc_age_coarse_cd4 = parse_data<int>(data, "hc_age_coarse_cd4", opts.p_idx_hiv_first_adult),
      .adult_female_infections = parse_data<real_type>(data, "adult_female_infections", opts.p_fertility_age_groups, proj_years),
      .adult_female_hivnpop = parse_data<real_type>(data, "hivnpop", opts.p_fertility_age_groups, proj_years),
      .total_births = parse_data<real_type>(data, "total_births", proj_years),
      .ctx_effect = parse_data<real_type>(data, "ctx_effect", 3),
      .hc_art_start = Rcpp::as<real_type>(data["hc_art_start"]),
      .local_adj_factor = Rcpp::as<real_type>(data["laf"])
    };
  };

  static constexpr int output_count = 13;

  static int build_output(
    Rcpp::List& ret,
    Rcpp::CharacterVector& names,
    int index,
    const Config::OutputState& state,
    const size_t& output_years
  ) {
    Rcpp::NumericVector r_hc1_hiv_pop(SS::hc1DS * SS::hcTT * SS::hc1AG * SS::NS * output_years);
    r_hc1_hiv_pop.attr("dim") = Rcpp::IntegerVector::create(SS::hc1DS, SS::hcTT, SS::hc1AG, SS::NS, output_years);
    std::copy_n(state.hc1_hiv_pop.data(), state.hc1_hiv_pop.size(), REAL(r_hc1_hiv_pop));
    names[index + 0] = "hc1_hiv_pop";
    ret[index + 0] = r_hc1_hiv_pop;
    Rcpp::NumericVector r_hc2_hiv_pop(SS::hc2DS * SS::hcTT * SS::hc2AG * SS::NS * output_years);
    r_hc2_hiv_pop.attr("dim") = Rcpp::IntegerVector::create(SS::hc2DS, SS::hcTT, SS::hc2AG, SS::NS, output_years);
    std::copy_n(state.hc2_hiv_pop.data(), state.hc2_hiv_pop.size(), REAL(r_hc2_hiv_pop));
    names[index + 1] = "hc2_hiv_pop";
    ret[index + 1] = r_hc2_hiv_pop;
    Rcpp::NumericVector r_hc1_art_pop(SS::hTS * SS::hc1DS * SS::hc1AG * SS::NS * output_years);
    r_hc1_art_pop.attr("dim") = Rcpp::IntegerVector::create(SS::hTS, SS::hc1DS, SS::hc1AG, SS::NS, output_years);
    std::copy_n(state.hc1_art_pop.data(), state.hc1_art_pop.size(), REAL(r_hc1_art_pop));
    names[index + 2] = "hc1_art_pop";
    ret[index + 2] = r_hc1_art_pop;
    Rcpp::NumericVector r_hc2_art_pop(SS::hTS * SS::hc2DS * SS::hc2AG * SS::NS * output_years);
    r_hc2_art_pop.attr("dim") = Rcpp::IntegerVector::create(SS::hTS, SS::hc2DS, SS::hc2AG, SS::NS, output_years);
    std::copy_n(state.hc2_art_pop.data(), state.hc2_art_pop.size(), REAL(r_hc2_art_pop));
    names[index + 3] = "hc2_art_pop";
    ret[index + 3] = r_hc2_art_pop;
    Rcpp::NumericVector r_hc1_noart_aids_deaths(SS::hc1DS * SS::hcTT * SS::hc1AG * SS::NS * output_years);
    r_hc1_noart_aids_deaths.attr("dim") = Rcpp::IntegerVector::create(SS::hc1DS, SS::hcTT, SS::hc1AG, SS::NS, output_years);
    std::copy_n(state.hc1_noart_aids_deaths.data(), state.hc1_noart_aids_deaths.size(), REAL(r_hc1_noart_aids_deaths));
    names[index + 4] = "hc1_noart_aids_deaths";
    ret[index + 4] = r_hc1_noart_aids_deaths;
    Rcpp::NumericVector r_hc2_noart_aids_deaths(SS::hc2DS * SS::hcTT * SS::hc2AG * SS::NS * output_years);
    r_hc2_noart_aids_deaths.attr("dim") = Rcpp::IntegerVector::create(SS::hc2DS, SS::hcTT, SS::hc2AG, SS::NS, output_years);
    std::copy_n(state.hc2_noart_aids_deaths.data(), state.hc2_noart_aids_deaths.size(), REAL(r_hc2_noart_aids_deaths));
    names[index + 5] = "hc2_noart_aids_deaths";
    ret[index + 5] = r_hc2_noart_aids_deaths;
    Rcpp::NumericVector r_hc1_art_aids_deaths(SS::hTS * SS::hc1DS * SS::hc1AG * SS::NS * output_years);
    r_hc1_art_aids_deaths.attr("dim") = Rcpp::IntegerVector::create(SS::hTS, SS::hc1DS, SS::hc1AG, SS::NS, output_years);
    std::copy_n(state.hc1_art_aids_deaths.data(), state.hc1_art_aids_deaths.size(), REAL(r_hc1_art_aids_deaths));
    names[index + 6] = "hc1_art_aids_deaths";
    ret[index + 6] = r_hc1_art_aids_deaths;
    Rcpp::NumericVector r_hc2_art_aids_deaths(SS::hTS * SS::hc2DS * SS::hc2AG * SS::NS * output_years);
    r_hc2_art_aids_deaths.attr("dim") = Rcpp::IntegerVector::create(SS::hTS, SS::hc2DS, SS::hc2AG, SS::NS, output_years);
    std::copy_n(state.hc2_art_aids_deaths.data(), state.hc2_art_aids_deaths.size(), REAL(r_hc2_art_aids_deaths));
    names[index + 7] = "hc2_art_aids_deaths";
    ret[index + 7] = r_hc2_art_aids_deaths;
    Rcpp::NumericVector r_hc_art_init(SS::hcAG_coarse * output_years);
    r_hc_art_init.attr("dim") = Rcpp::IntegerVector::create(SS::hcAG_coarse, output_years);
    std::copy_n(state.hc_art_init.data(), state.hc_art_init.size(), REAL(r_hc_art_init));
    names[index + 8] = "hc_art_init";
    ret[index + 8] = r_hc_art_init;
    Rcpp::NumericVector r_hc_art_need_init(SS::hc1DS * SS::hcTT * SS::hcAG_end * SS::NS * output_years);
    r_hc_art_need_init.attr("dim") = Rcpp::IntegerVector::create(SS::hc1DS, SS::hcTT, SS::hcAG_end, SS::NS, output_years);
    std::copy_n(state.hc_art_need_init.data(), state.hc_art_need_init.size(), REAL(r_hc_art_need_init));
    names[index + 9] = "hc_art_need_init";
    ret[index + 9] = r_hc_art_need_init;
    Rcpp::NumericVector r_hiv_births(output_years);
    r_hiv_births.attr("dim") = Rcpp::IntegerVector::create(output_years);
    std::copy_n(state.hiv_births.data(), state.hiv_births.size(), REAL(r_hiv_births));
    names[index + 10] = "hiv_births";
    ret[index + 10] = r_hiv_births;
    Rcpp::NumericVector r_ctx_need(output_years);
    r_ctx_need.attr("dim") = Rcpp::IntegerVector::create(output_years);
    std::copy_n(state.ctx_need.data(), state.ctx_need.size(), REAL(r_ctx_need));
    names[index + 11] = "ctx_need";
    ret[index + 11] = r_ctx_need;
    Rcpp::NumericVector r_infection_by_type(SS::hcTT * SS::hc1AG * SS::NS * output_years);
    r_infection_by_type.attr("dim") = Rcpp::IntegerVector::create(SS::hcTT, SS::hc1AG, SS::NS, output_years);
    std::copy_n(state.infection_by_type.data(), state.infection_by_type.size(), REAL(r_infection_by_type));
    names[index + 12] = "infection_by_type";
    ret[index + 12] = r_infection_by_type;
    return index + output_count;
  };
};


}
