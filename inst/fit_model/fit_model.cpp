#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <vector>

#include <frogger.hpp>
#include <serialize_eigen.hpp>

template<typename T, int rank>
Eigen::TensorMap <Eigen::Tensor<T, rank>> tensor_to_tensor_map(Eigen::Tensor <T, rank> &d) {
  return Eigen::TensorMap < Eigen::Tensor < T, rank >> (d.data(), d.dimensions());
}

int main(int argc, char *argv[]) {
  if (argc < 4) {
    std::cout <<
              "Usage: fit_model <sim_years> <hiv_steps_per_year> <intput_dir> <output_dir>" <<
              std::endl;
    return 1;
  }

  int sim_years = atoi(argv[1]);
  int hiv_steps_per_year = atoi(argv[2]);
  std::string input_dir = argv[3];
  std::string output_dir = argv[4];

  std::string input_abs = std::filesystem::absolute(input_dir);
  if (!std::filesystem::exists(input_abs)) {
    std::cout << "Input dir '" << input_dir << "' does not exist." << std::endl;
    return 1;
  }

  std::string output_abs = std::filesystem::absolute(output_dir);
  if (!std::filesystem::exists(output_abs)) {
    if (std::filesystem::create_directory(output_abs)) {
      std::cout << "Created output directory '" << output_abs << "'" << std::endl;
    } else {
      std::cout << "Failed to create output directory '" << output_abs << "'" << std::endl;
    }
  } else {
    std::cout << "Writing to existing output directory " << output_abs << "'" << std::endl;
  }

  if (sim_years > 60) {
    std::cout << "Running to max no of sim years: 60\n" << std::endl;
    sim_years = 60;
  }
  if (hiv_steps_per_year > 10) {
    std::cout << "Running max no of HIV steps per years: 10" << std::endl;
    hiv_steps_per_year = 10;
  }

  const double dt = (1.0 / hiv_steps_per_year);
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
  const int time_art_start = 30;

  // Set working dir to read files
  std::filesystem::path old_dir = std::filesystem::current_path();
  std::filesystem::current_path(input_abs);

  // Only fine-grained ages at first
  leapfrog::Tensor1<int> age_groups_hiv_span_data = serialize::deserialize_tensor<int, 1>(
      std::string("hAG_SPAN_full"));
  const leapfrog::TensorMap1<int> age_groups_hiv_span = tensor_to_tensor_map<int, 1>(age_groups_hiv_span_data);
  int age_groups_hiv = static_cast<int>(age_groups_hiv_span.size());
  int age_groups_hiv_15plus = age_groups_hiv;
  const int scale_cd4_mortality = 1;
  int hIDX_15PLUS = 0;
  const double art_alloc_mxweight = 0.2;

  leapfrog::Tensor1<int> v = serialize::deserialize_tensor<int, 1>(std::string("artcd4elig_idx"));
  for (int i = 0; i <= sim_years; ++i) {
    // 0-based indexing in C++ vs 1-based indexing in R
    v(i) = v[i] - 1;
  }
  const leapfrog::TensorMap1<int> artcd4elig_idx = tensor_to_tensor_map<int, 1>(v);

  leapfrog::Tensor1<double> h(treatment_stages - 1);
  for (int i = 0; i < treatment_stages - 1; ++i) {
    h(i) = 0.5;
  }
  const leapfrog::TensorMap1<double> h_art_stage_dur = tensor_to_tensor_map<double, 1>(h);

  leapfrog::Tensor2<double> base_pop_data = serialize::deserialize_tensor<double, 2>(
      std::string("basepop"));
  const leapfrog::TensorMap2<double> base_pop = tensor_to_tensor_map<double, 2>(base_pop_data);
  leapfrog::Tensor3<double> survival_data = serialize::deserialize_tensor<double, 3>(
      std::string("survival"));
  const leapfrog::TensorMap3<double> survival = tensor_to_tensor_map<double, 3>(survival_data);
  leapfrog::Tensor3<double> net_migration_data = serialize::deserialize_tensor<double, 3>(
      std::string("net_migration"));
  const leapfrog::TensorMap3<double> net_migration = tensor_to_tensor_map<double, 3>(net_migration_data);
  leapfrog::Tensor2<double> age_sex_fertility_ratio_data = serialize::deserialize_tensor<double, 2>(
      std::string("age_sex_fertility_ratio"));
  const leapfrog::TensorMap2<double> age_sex_fertility_ratio = tensor_to_tensor_map<double, 2>(
      age_sex_fertility_ratio_data);
  leapfrog::Tensor2<double> births_sex_prop_data = serialize::deserialize_tensor<double, 2>(
      std::string("births_sex_prop"));
  const leapfrog::TensorMap2<double> births_sex_prop = tensor_to_tensor_map<double, 2>(births_sex_prop_data);
  leapfrog::Tensor1<double> incidence_rate_data = serialize::deserialize_tensor<double, 1>(
      std::string("incidence_rate"));
  const leapfrog::TensorMap1<double> incidence_rate = tensor_to_tensor_map<double, 1>(incidence_rate_data);
  leapfrog::Tensor3<double> incidence_relative_risk_age_data = serialize::deserialize_tensor<double, 3>(
      std::string("incidence_rate_relative_risk_age"));
  const leapfrog::TensorMap3<double> incidence_relative_risk_age = tensor_to_tensor_map<double, 3>(
      incidence_relative_risk_age_data);
  leapfrog::Tensor1<double> incidence_relative_risk_sex_data = serialize::deserialize_tensor<double, 1>(
      std::string("incidence_rate_relative_risk_sex"));
  const leapfrog::TensorMap1<double> incidence_relative_risk_sex = tensor_to_tensor_map<double, 1>(
      incidence_relative_risk_sex_data);
  leapfrog::Tensor3<double> cd4_mortality_data = serialize::deserialize_tensor<double, 3>(
      std::string("cd4_mortality_full"));
  const leapfrog::TensorMap3<double> cd4_mortality = tensor_to_tensor_map<double, 3>(cd4_mortality_data);
  leapfrog::Tensor3<double> cd4_progression_data = serialize::deserialize_tensor<double, 3>(
      std::string("cd4_progression_full"));
  const leapfrog::TensorMap3<double> cd4_progression = tensor_to_tensor_map<double, 3>(cd4_progression_data);
  leapfrog::Tensor3<double> cd4_initdist_data = serialize::deserialize_tensor<double, 3>(
      std::string("cd4_initdist_full"));
  const leapfrog::TensorMap3<double> cd4_initdist = tensor_to_tensor_map<double, 3>(cd4_initdist_data);
  const leapfrog::TensorMap1<int> hiv_age_groups_span = tensor_to_tensor_map<int, 1>(age_groups_hiv_span_data);
  leapfrog::Tensor4<double> art_mortality_data = serialize::deserialize_tensor<double, 4>(
      std::string("art_mortality_full"));
  const leapfrog::TensorMap4<double> art_mortality = tensor_to_tensor_map<double, 4>(art_mortality_data);
  leapfrog::Tensor2<double> artmx_timerr_data = serialize::deserialize_tensor<double, 2>(
      std::string("artmx_timerr"));
  const leapfrog::TensorMap2<double> artmx_timerr = tensor_to_tensor_map<double, 2>(artmx_timerr_data);
  leapfrog::Tensor1<double> art_dropout_data = serialize::deserialize_tensor<double, 1>(
      std::string("art_dropout"));
  const leapfrog::TensorMap1<double> art_dropout = tensor_to_tensor_map<double, 1>(art_dropout_data);
  Eigen::Tensor<double, 2> art15plus_num_data = serialize::deserialize_tensor<double, 2>(
      std::string("art15plus_num"));
  const leapfrog::TensorMap2<double> art15plus_num = tensor_to_tensor_map<double, 2>(art15plus_num_data);
  leapfrog::Tensor2<int> art15plus_isperc_data = serialize::deserialize_tensor<int, 2>(
      std::string("art15plus_isperc"));
  const leapfrog::TensorMap2<int> art15plus_isperc = tensor_to_tensor_map<int, 2>(art15plus_isperc_data);

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
                                               hiv_steps_per_year,
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

  auto output_state = leapfrog::run_model(sim_years, params);
  std::cout << "Fit complete" << std::endl;

  std::filesystem::path out_path(output_abs);
  std::filesystem::path total_population_path = out_path / "total_population";
  serialize::serialize_tensor<double, 2>(output_state.total_population, total_population_path);
  std::filesystem::path births_path = out_path / "births";

  // Births is just a double so write it out
  std::ofstream dest(births_path);
  dest << output_state.births << std::endl;
  dest.close();

  std::filesystem::path natural_deaths_path = out_path / "natural_deaths";
  serialize::serialize_tensor<double, 2>(output_state.natural_deaths, natural_deaths_path);
  std::filesystem::path hiv_population_path = out_path / "hiv_population";
  serialize::serialize_tensor<double, 2>(output_state.hiv_population, hiv_population_path);
  std::filesystem::path hiv_natural_deaths_path = out_path / "hiv_natural_deaths";
  serialize::serialize_tensor<double, 2>(output_state.hiv_natural_deaths, hiv_natural_deaths_path);
  std::filesystem::path hiv_strat_adult_path = out_path / "hiv_strat_adult";
  serialize::serialize_tensor<double, 3>(output_state.hiv_strat_adult, hiv_strat_adult_path);
  std::filesystem::path art_strat_adult_path = out_path / "art_strat_adult";
  serialize::serialize_tensor<double, 4>(output_state.art_strat_adult, art_strat_adult_path);
  std::filesystem::path aids_deaths_no_art_path = out_path / "aids_deaths_no_art";
  serialize::serialize_tensor<double, 3>(output_state.aids_deaths_no_art, aids_deaths_no_art_path);
  std::filesystem::path infections_path = out_path / "infections";
  serialize::serialize_tensor<double, 2>(output_state.infections, infections_path);
  std::filesystem::path aids_deaths_art_path = out_path / "aids_deaths_art";
  serialize::serialize_tensor<double, 4>(output_state.aids_deaths_art, aids_deaths_art_path);
  std::filesystem::path art_initiation_path = out_path / "art_initiation";
  serialize::serialize_tensor<double, 3>(output_state.art_initiation, art_initiation_path);
  std::filesystem::path hiv_deaths_path = out_path / "hiv_deaths";
  serialize::serialize_tensor<double, 2>(output_state.hiv_deaths, hiv_deaths_path);

  return 0;
}
