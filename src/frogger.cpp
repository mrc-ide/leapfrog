#include <Rcpp.h>

#include "frogger.hpp"
#include "types.hpp"

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

int get_age_groups_hiv(const std::string hiv_age_stratification) {
  int age_groups_hiv;
  if (hiv_age_stratification == "full") {
    age_groups_hiv = 66;
  } else {
    // We've already validated that hiv_age_stratification is one of "full" or "coarse"
    age_groups_hiv = 9;
  }
  return age_groups_hiv;
}

template<typename... Args>
auto parse_data_int(const Rcpp::List data, const std::string& key, Args... dims) {
  constexpr std::size_t rank = sizeof...(dims);
  Eigen::array<int, rank> dimensions{ static_cast<int>(dims)... };

  int product = 1;
  for (size_t i = 0; i < rank; ++i) {
    product *= dimensions[i];
  }
  SEXP array_data = data[key];
  // In cases where the input data has project years we might not use all of it model fit
  // So we can take create a Map over a smaller slice of the data
  // As long as this is true we can be confident we're not referencing invalid memory
  if (LENGTH(array_data) < product) {
    Rcpp::stop("Invalid size of data for '%s', expected %d got %d",
               key,
               product,
               LENGTH(array_data));
  }

  return Eigen::TensorMap<Eigen::Tensor<int, rank>>(INTEGER(array_data), static_cast<int>(dims)...);
}

template<typename... Args>
auto parse_data_double(const Rcpp::List data, const std::string& key, Args... dims) {
  constexpr std::size_t rank = sizeof...(dims);
  Eigen::array<int, rank> dimensions{ static_cast<int>(dims)... };

  int product = 1;
  for (size_t i = 0; i < rank; ++i) {
    product *= dimensions[i];
  }
  SEXP array_data = data[key];
  // In cases where the input data has project years we might not use all of it model fit
  // So we can take create a Map over a smaller slice of the data
  // As long as this is true we can be confident we're not referencing invalid memory
  if (LENGTH(array_data) < product) {
    Rcpp::stop("Invalid size of data for '%s', expected %d got %d.",
               key,
               product,
               LENGTH(array_data));
  }

  return Eigen::TensorMap<Eigen::Tensor<double, rank>>(REAL(array_data), static_cast<int>(dims)...);
}

std::vector<int> parse_output_steps(Rcpp::NumericVector output_steps, int proj_years) {
  return Rcpp::as<std::vector<int>>(output_steps);
}

void validate_stratification(const std::string stratification) {
  if (!(stratification == "full" || stratification == "coarse")) {
    Rcpp::stop(
        "Invalid HIV age stratification must be 'full' or 'coarse' got '%s'.",
        stratification);
  }
}

// [[Rcpp::export]]
Rcpp::List run_base_model(const Rcpp::List data,
                          const Rcpp::List projection_parameters,
                          SEXP sim_years,
                          SEXP hiv_steps_per_year,
                          Rcpp::NumericVector output_steps,
                          std::string hiv_age_stratification = "full") {
  validate_stratification(hiv_age_stratification);
  const int proj_years = get_simulation_years(data, sim_years);
  const std::vector<int> save_steps = parse_output_steps(output_steps, proj_years);
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
  const int age_groups_hiv = get_age_groups_hiv(hiv_age_stratification);

  int age_groups_hiv_15plus = age_groups_hiv;
  const int scale_cd4_mortality =
      Rcpp::as<int>(projection_parameters["scale_cd4_mort"]);
  int hIDX_15PLUS = 0;
  const double art_alloc_mxweight = Rcpp::as<double>(projection_parameters["art_alloc_mxweight"]);

  std::function<std::string(const std::string&)> build_prop_name = [hiv_age_stratification](const std::string& property) {
    return property + "_" + hiv_age_stratification;
  };

  const leapfrog::TensorMap1<int> age_groups_hiv_span = parse_data_int(
      projection_parameters, build_prop_name("hAG_SPAN"), age_groups_hiv);
  const leapfrog::TensorMap2<double> base_pop = parse_data_double(data, "basepop", age_groups_pop, num_genders);

  // Survival has size age_groups_pop + 1 as this is the probability of
  // surviving between ages, so from 0 to 1, 1 to 2, ..., 79 to 80+ and
  // 80+ to 80+
  const leapfrog::TensorMap3<double> survival = parse_data_double(data, "Sx",
                                                                  age_groups_pop + 1, num_genders, proj_years);
  const leapfrog::TensorMap3<double> net_migration = parse_data_double(data, "netmigr_adj", age_groups_pop,
                                                                       num_genders, proj_years);
  const leapfrog::TensorMap2<double> age_sex_fertility_ratio = parse_data_double(data, "asfr",
                                                                                 age_groups_fert, proj_years);
  const leapfrog::TensorMap2<double> births_sex_prop = parse_data_double(data, "births_sex_prop",
                                                                         num_genders, proj_years);
  const leapfrog::TensorMap1<double> incidence_rate = parse_data_double(projection_parameters, "incidinput",
                                                                        proj_years);
  const leapfrog::TensorMap3<double> incidence_relative_risk_age = parse_data_double(
      projection_parameters, "incrr_age",age_groups_pop - hiv_adult_first_age_group, num_genders, proj_years);
  const leapfrog::TensorMap1<double> incidence_relative_risk_sex = parse_data_double(projection_parameters, "incrr_sex",
                                                                                     proj_years);
  const leapfrog::TensorMap3<double> cd4_mortality = parse_data_double(
      projection_parameters, build_prop_name("cd4_mort"), disease_stages, age_groups_hiv, num_genders);
  const leapfrog::TensorMap3<double> cd4_progression = parse_data_double(
      projection_parameters, build_prop_name("cd4_prog"), disease_stages - 1, age_groups_hiv, num_genders);

  leapfrog::Tensor1<int> artcd4elig_idx = parse_data_int(projection_parameters, "artcd4elig_idx", proj_years + 1);
  for (int i = 0; i <= proj_years; ++i) {
    // 0-based indexing in C++ vs 1-based indexing in R
    artcd4elig_idx(i) = artcd4elig_idx(i) - 1;
  }

  const leapfrog::TensorMap3<double> cd4_initdist = parse_data_double(
      projection_parameters, build_prop_name("cd4_initdist"), disease_stages, age_groups_hiv, num_genders);
  const leapfrog::TensorMap1<int> hiv_age_groups_span = parse_data_int(
      projection_parameters, build_prop_name("hAG_SPAN"), age_groups_hiv);
  const leapfrog::TensorMap4<double> art_mortality = parse_data_double(
      projection_parameters, build_prop_name("art_mort"), treatment_stages, disease_stages, age_groups_hiv,
      num_genders);
  const leapfrog::TensorMap2<double> artmx_timerr = parse_data_double(projection_parameters, "artmx_timerr",
                                                                      treatment_stages, proj_years);
  leapfrog::Tensor1<double> h_art_stage_dur(treatment_stages - 1);
  for (int i = 0; i < treatment_stages - 1; ++i) {
    h_art_stage_dur(i) = 0.5;
  }
  const leapfrog::TensorMap1<double> art_dropout = parse_data_double(projection_parameters, "art_dropout", proj_years);
  const leapfrog::TensorMap2<double> art15plus_num = parse_data_double(projection_parameters, "art15plus_num",
                                                                       num_genders, proj_years);
  const leapfrog::TensorMap2<int> art15plus_isperc = parse_data_int(projection_parameters, "art15plus_isperc",
                                                                    num_genders, proj_years);

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

  auto state = leapfrog::run_model(proj_years, save_steps, params);

  size_t output_years = save_steps.size();
  Rcpp::NumericVector r_total_population(age_groups_pop * num_genders * output_years);
  Rcpp::NumericVector r_births(output_years);
  Rcpp::NumericVector r_natural_deaths(age_groups_pop * num_genders * output_years);
  Rcpp::NumericVector r_hiv_population(age_groups_pop * num_genders * output_years);
  Rcpp::NumericVector r_hiv_natural_deaths(age_groups_pop * num_genders * output_years);
  Rcpp::NumericVector r_hiv_strat_adult(disease_stages * age_groups_hiv *
                                        num_genders * output_years);
  Rcpp::NumericVector r_art_strat_adult(treatment_stages * disease_stages *
                                        age_groups_hiv * num_genders * output_years);
  Rcpp::NumericVector r_aids_deaths_no_art(disease_stages * age_groups_hiv * num_genders * output_years);
  Rcpp::NumericVector r_infections(age_groups_pop * num_genders * output_years);
  Rcpp::NumericVector r_aids_deaths_art(treatment_stages * disease_stages * age_groups_hiv *
                                        num_genders * output_years);
  Rcpp::NumericVector r_art_initiation(disease_stages * age_groups_hiv * num_genders * output_years);
  Rcpp::NumericVector r_hiv_deaths(age_groups_pop * num_genders * output_years);

  r_total_population.attr("dim") =
      Rcpp::NumericVector::create(age_groups_pop, num_genders, output_years);
  r_natural_deaths.attr("dim") =
      Rcpp::NumericVector::create(age_groups_pop, num_genders, output_years);
  r_hiv_population.attr("dim") =
      Rcpp::NumericVector::create(age_groups_pop, num_genders, output_years);
  r_hiv_natural_deaths.attr("dim") =
      Rcpp::NumericVector::create(age_groups_pop, num_genders, output_years);
  r_hiv_strat_adult.attr("dim") =
      Rcpp::NumericVector::create(disease_stages, age_groups_hiv, num_genders, output_years);
  r_art_strat_adult.attr("dim") = Rcpp::NumericVector::create(
      treatment_stages, disease_stages, age_groups_hiv, num_genders, output_years);
  r_births.attr("dim") = Rcpp::NumericVector::create(output_years);
  r_aids_deaths_no_art.attr("dim") = Rcpp::NumericVector::create(
      disease_stages, age_groups_hiv, num_genders, output_years);
  r_infections.attr("dim") = Rcpp::NumericVector::create(age_groups_pop, num_genders, output_years);
  r_aids_deaths_art.attr("dim") = Rcpp::NumericVector::create(
      treatment_stages, disease_stages, age_groups_hiv, num_genders, output_years);
  r_art_initiation.attr("dim") = Rcpp::NumericVector::create(
      disease_stages, age_groups_hiv, num_genders, output_years);
  r_hiv_deaths.attr("dim") = Rcpp::NumericVector::create(age_groups_pop, num_genders, output_years);

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
  std::copy_n(state.births.data(), state.births.size(),
              REAL(r_births));
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
