#include <frogger.hpp>
#include <vector>
#include <fstream>

leapfrog::Parameters<double> build_parameters(int sim_years, int hiv_steps_per_year) {
  if (sim_years > 60) {
    std::cout <<  "Running to max no of sim years: 60"
    sim_years = 60;
  }

  if (hiv_steps_per_year > 10) {
    std::cout <<  "Running max no of HIV steps per years: 10"
    hiv_steps_per_year = 10;
  }
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
}

int main(int argc, char *argv[]) {

  auto params = build_parameters(61, 10);
  auto state = leapfrog::run_model(61, params);

  return 0;
}
