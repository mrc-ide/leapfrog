#include <Rcpp.h>

#include "frogger.hpp"
#include "types.hpp"
#include "serialize_eigen.hpp"

int get_simulation_years(const Rcpp::List demp, SEXP r_sim_years) {
  Rcpp::NumericVector Sx = demp["Sx"];
  Rcpp::Dimension d = Sx.attr("dim");
  // Simulation initialises state from first years input data (index 0)
  // then runs for each year simulating this years (i) data using previous years
  // state (i - 1) and this years input data (i). So -1 off index for max years
  // to simulate as index 0 used for initial state
  const int max_sim_years = d[2] - 1;
  if (r_sim_years == R_NilValue) {
    return max_sim_years;
  }
  auto sim_years = INTEGER(r_sim_years)[0];
  if (sim_years > max_sim_years) {
    Rcpp::stop("No of years > max years of " + std::to_string(max_sim_years));
  }
  return sim_years;
}

int get_hiv_steps_per_year(SEXP r_hiv_steps_per_year) {
  int hiv_steps_per_year;
  if (r_hiv_steps_per_year == R_NilValue) {
    hiv_steps_per_year = 10;
  } else {
    hiv_steps_per_year = INTEGER(r_hiv_steps_per_year)[0];
  }
  return hiv_steps_per_year;
}

leapfrog::TensorMap1<int> get_age_groups_hiv_span(const Rcpp::List projection_parameters,
                                                  std::string hiv_age_stratification) {
  int age_groups_hiv;
  SEXP data;
  if (hiv_age_stratification == "full") {
    age_groups_hiv = 66;
    data = projection_parameters["hAG_SPAN_full"];
  } else if (hiv_age_stratification == "coarse") {
    age_groups_hiv = 9;
    data = projection_parameters["hAG_SPAN_coarse"];
  } else {
    Rcpp::stop(
        "Invalid HIV age stratification must be 'full' or 'coarse' got '%s'.",
        hiv_age_stratification);
  }
  if (LENGTH(data) != age_groups_hiv) {
    Rcpp::stop("Invalid size of data, expected %d got %d", age_groups_hiv,
               LENGTH(data));
  }
  leapfrog::TensorMap1<int> age_groups_hiv_span(INTEGER(data), age_groups_hiv);
  return age_groups_hiv_span;
}

// [[Rcpp::export]]
Rcpp::List run_base_model(const Rcpp::List data,
                          const Rcpp::List projection_parameters,
                          SEXP sim_years,
                          SEXP hiv_steps_per_year,
                          std::string hiv_age_stratification = "full") {
  const int proj_years = get_simulation_years(data, sim_years);
  const int hiv_steps = get_hiv_steps_per_year(hiv_steps_per_year);
  const double dt = (1.0 / hiv_steps);
  const int num_genders = 2;
  const int age_groups_pop = 81;
  const int fertility_first_age_group = 15;
  const int age_groups_fert = 35;
  const int disease_stages = 7;
  const int treatment_stages = 3;
  const int hiv_adult_first_age_group = 15;
  const int adult_incidence_first_age_group = hiv_adult_first_age_group;
  // Hardcoded 15-49 for now (35 groups within this band)
  const int pAG_INCIDPOP = 35;
  // 0-based indexing vs R 1-based
  const int time_art_start =
      Rcpp::as<int>(projection_parameters["t_ART_start"]) - 1;
  const leapfrog::TensorMap1<int> age_groups_hiv_span =
      get_age_groups_hiv_span(projection_parameters, hiv_age_stratification);
  int age_groups_hiv = static_cast<int>(age_groups_hiv_span.size());
  int age_groups_hiv_15plus = age_groups_hiv;
  const int scale_cd4_mortality =
      Rcpp::as<int>(projection_parameters["scale_cd4_mort"]);
  int hIDX_15PLUS = 0;
  const double art_alloc_mxweight = Rcpp::as<double>(projection_parameters["art_alloc_mxweight"]);

  const leapfrog::TensorMap2<double> base_pop(REAL(data["basepop"]), age_groups_pop,
                                              num_genders);
  // Survival has size age_groups_pop + 1 as this is the probability of
  // surviving between ages, so from 0 to 1, 1 to 2, ..., 79 to 80+ and
  // 80+ to 80+
  const leapfrog::TensorMap3<double> survival(REAL(data["Sx"]), age_groups_pop + 1, num_genders,
                                              proj_years);
  const leapfrog::TensorMap3<double> net_migration(REAL(data["netmigr_adj"]), age_groups_pop,
                                                   num_genders, proj_years);
  const leapfrog::TensorMap2<double> age_sex_fertility_ratio(REAL(data["asfr"]),
                                                             age_groups_fert, proj_years);
  const leapfrog::TensorMap2<double> births_sex_prop(REAL(data["births_sex_prop"]), num_genders,
                                                     proj_years);
  const leapfrog::TensorMap1<double> incidence_rate(REAL(projection_parameters["incidinput"]), proj_years);
  const leapfrog::TensorMap3<double> incidence_relative_risk_age(REAL(projection_parameters["incrr_age"]),
                                                                 age_groups_pop - hiv_adult_first_age_group,
                                                                 num_genders, proj_years);
  const leapfrog::TensorMap1<double> incidence_relative_risk_sex(REAL(projection_parameters["incrr_sex"]), proj_years);
  const leapfrog::TensorMap3<double> cd4_mortality(REAL(projection_parameters["cd4_mort_full"]),
                                                   disease_stages, age_groups_hiv, num_genders);
  const leapfrog::TensorMap3<double> cd4_progression(REAL(projection_parameters["cd4_prog_full"]),
                                                     disease_stages - 1, age_groups_hiv, num_genders);
  Rcpp::IntegerVector v = Rcpp::as<Rcpp::IntegerVector>(projection_parameters["artcd4elig_idx"]);
  leapfrog::Tensor1<int> artcd4elig_idx(proj_years + 1);
  for (int i = 0; i <= proj_years; ++i) {
    // 0-based indexing in C++ vs 1-based indexing in R
    artcd4elig_idx(i) = v[i] - 1;
  }
  const leapfrog::TensorMap3<double> cd4_initdist(REAL(projection_parameters["cd4_initdist_full"]), disease_stages,
                                                  age_groups_hiv, num_genders);
  const leapfrog::TensorMap1<int> hiv_age_groups_span(INTEGER(projection_parameters["hAG_SPAN_full"]), age_groups_hiv);
  const leapfrog::TensorMap4<double> art_mortality(REAL(projection_parameters["art_mort_full"]), treatment_stages,
                                                   disease_stages, age_groups_hiv, num_genders);
  const leapfrog::TensorMap2<double> artmx_timerr(REAL(projection_parameters["artmx_timerr"]), treatment_stages,
                                                  proj_years);
  leapfrog::Tensor1<double> h_art_stage_dur(treatment_stages - 1);
  for (int i = 0; i < treatment_stages - 1; ++i) {
    h_art_stage_dur(i) = 0.5;
  }
  const leapfrog::TensorMap1<double> art_dropout(REAL(projection_parameters["art_dropout"]), proj_years);
  const leapfrog::TensorMap2<double> art15plus_num(REAL(projection_parameters["art15plus_num"]), num_genders,
                                                   proj_years);
  const leapfrog::TensorMap2<int> art15plus_isperc(INTEGER(projection_parameters["art15plus_isperc"]), num_genders,
                                                   proj_years);

  const leapfrog::Parameters<double> params = {num_genders,
                                               age_groups_pop,
                                               fertility_first_age_group,
                                               age_groups_fert,
                                               age_groups_hiv,
                                               age_groups_hiv_15plus,
                                               disease_stages,
                                               hiv_adult_first_age_group,
                                               treatment_stages,
                                               time_art_start,
                                               adult_incidence_first_age_group,
                                               pAG_INCIDPOP,
                                               hiv_steps,
                                               dt,
                                               scale_cd4_mortality,
                                               hIDX_15PLUS,
                                               art_alloc_mxweight,
                                               age_groups_hiv_span,
                                               incidence_rate,
                                               base_pop,
                                               survival,
                                               net_migration,
                                               age_sex_fertility_ratio,
                                               births_sex_prop,
                                               incidence_relative_risk_age,
                                               incidence_relative_risk_sex,
                                               cd4_mortality,
                                               cd4_progression,
                                               artcd4elig_idx,
                                               cd4_initdist,
                                               hiv_age_groups_span,
                                               art_mortality,
                                               artmx_timerr,
                                               h_art_stage_dur,
                                               art_dropout,
                                               art15plus_num,
                                               art15plus_isperc};

  auto state = leapfrog::run_model(proj_years, params);

  Rcpp::NumericVector r_total_population(age_groups_pop * num_genders);
  Rcpp::NumericVector r_births(1);
  Rcpp::NumericVector r_natural_deaths(age_groups_pop * num_genders);
  Rcpp::NumericVector r_hiv_population(age_groups_pop * num_genders);
  Rcpp::NumericVector r_hiv_natural_deaths(age_groups_pop * num_genders);
  Rcpp::NumericVector r_hiv_strat_adult(disease_stages * age_groups_hiv *
                                        num_genders);
  Rcpp::NumericVector r_art_strat_adult(treatment_stages * disease_stages *
                                        age_groups_hiv * num_genders);
  Rcpp::NumericVector r_aids_deaths_no_art(disease_stages * age_groups_hiv * num_genders);
  Rcpp::NumericVector r_infections(age_groups_pop * num_genders);
  Rcpp::NumericVector r_aids_deaths_art(treatment_stages * disease_stages * age_groups_hiv * num_genders);
  Rcpp::NumericVector r_art_initiation(disease_stages * age_groups_hiv * num_genders);
  Rcpp::NumericVector r_hiv_deaths(age_groups_pop * num_genders);

  r_total_population.attr("dim") =
      Rcpp::NumericVector::create(age_groups_pop, num_genders);
  r_natural_deaths.attr("dim") =
      Rcpp::NumericVector::create(age_groups_pop, num_genders);
  r_hiv_population.attr("dim") =
      Rcpp::NumericVector::create(age_groups_pop, num_genders);
  r_hiv_natural_deaths.attr("dim") =
      Rcpp::NumericVector::create(age_groups_pop, num_genders);
  r_hiv_strat_adult.attr("dim") =
      Rcpp::NumericVector::create(disease_stages, age_groups_hiv, num_genders);
  r_art_strat_adult.attr("dim") = Rcpp::NumericVector::create(
      treatment_stages, disease_stages, age_groups_hiv, num_genders);
  r_aids_deaths_no_art.attr("dim") = Rcpp::NumericVector::create(
      disease_stages, age_groups_hiv, num_genders);
  r_infections.attr("dim") = Rcpp::NumericVector::create(age_groups_pop, num_genders);
  r_aids_deaths_art.attr("dim") = Rcpp::NumericVector::create(
      treatment_stages, disease_stages, age_groups_hiv, num_genders);
  r_art_initiation.attr("dim") = Rcpp::NumericVector::create(
      disease_stages, age_groups_hiv, num_genders);
  r_hiv_deaths.attr("dim") = Rcpp::NumericVector::create(age_groups_pop, num_genders);

  std::copy_n(state.total_population.data(), state.total_population.size(),
              REAL(r_total_population));
  std::copy_n(state.natural_deaths.data(), state.natural_deaths.size(),
              REAL(r_natural_deaths));
  std::copy_n(state.hiv_population.data(), state.hiv_population.size(),
              REAL(r_hiv_population));
  std::copy_n(state.hiv_natural_deaths.data(), state.hiv_natural_deaths.size(),
              REAL(r_hiv_natural_deaths));
  std::copy_n(state.hiv_strat_adult.data(), state.hiv_strat_adult.size(),
              REAL(r_hiv_strat_adult));
  std::copy_n(state.art_strat_adult.data(), state.art_strat_adult.size(),
              REAL(r_art_strat_adult));
  std::copy_n(state.aids_deaths_no_art.data(), state.aids_deaths_no_art.size(),
              REAL(r_aids_deaths_no_art));
  std::copy_n(state.infections.data(), state.infections.size(),
              REAL(r_infections));
  std::copy_n(state.aids_deaths_art.data(), state.aids_deaths_art.size(),
              REAL(r_aids_deaths_art));
  std::copy_n(state.art_initiation.data(), state.art_initiation.size(),
              REAL(r_art_initiation));
  std::copy_n(state.hiv_deaths.data(), state.hiv_deaths.size(),
              REAL(r_hiv_deaths));
  REAL(r_births)[0] = state.births;

  Rcpp::List ret =
      Rcpp::List::create(Rcpp::_["total_population"] = r_total_population,
                         Rcpp::_["births"] = r_births,
                         Rcpp::_["natural_deaths"] = r_natural_deaths,
                         Rcpp::_["hiv_population"] = r_hiv_population,
                         Rcpp::_["hiv_natural_deaths"] = r_hiv_natural_deaths,
                         Rcpp::_["hiv_strat_adult"] = r_hiv_strat_adult,
                         Rcpp::_["art_strat_adult"] = r_art_strat_adult,
                         Rcpp::_["aids_deaths_no_art"] = r_aids_deaths_no_art,
                         Rcpp::_["infections"] = r_infections,
                         Rcpp::_["aids_deaths_art"] = r_aids_deaths_art,
                         Rcpp::_["art_initiation"] = r_art_initiation,
                         Rcpp::_["hiv_deaths"] = r_hiv_deaths);
  return ret;
}


// [[Rcpp::export]]
Rcpp::List serialize_vector(const Rcpp::List data, const std::string path1, const std::string path2) {
  const leapfrog::TensorMap2<double> test_2d(REAL(data["test_2d_double"]), 3, 4);
  const leapfrog::TensorMap3<int> test_3d(INTEGER(data["test_3d_int"]), 4, 3, 2);

  serialize::serialize_tensor_map(test_2d, path1);
  serialize::serialize_tensor_map(test_3d, path2);

  Rcpp::List ret =
      Rcpp::List::create(Rcpp::_["test_2d_double_path"] = path1,
                         Rcpp::_["test_3d_int_path"] = path2);
  return ret;
}

// [[Rcpp::export]]
Rcpp::List deserialize_vector(const std::string path1, const std::string path2) {

  const Eigen::Tensor<double, 2> test_2d = serialize::deserialize_tensor<double, 2>(path1);
  const Eigen::Tensor<int, 3> test_3d = serialize::deserialize_tensor<int, 3>(path2);

  Rcpp::NumericVector r_test_2d(3 * 4);
  Rcpp::NumericVector r_test_3d(4 * 3 * 2);

  r_test_2d.attr("dim") =
      Rcpp::NumericVector::create(3, 4);
  r_test_3d.attr("dim") =
      Rcpp::NumericVector::create(4, 3, 2);

  std::copy_n(test_2d.data(), test_2d.size(), REAL(r_test_2d));
  std::copy_n(test_3d.data(), test_3d.size(), INTEGER(r_test_3d));

  Rcpp::List ret =
      Rcpp::List::create(Rcpp::_["test_2d_double"] = r_test_2d,
                         Rcpp::_["test_3d_int"] = r_test_3d);
  return ret;
}
