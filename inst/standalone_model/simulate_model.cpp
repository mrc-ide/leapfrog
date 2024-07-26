#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <vector>

#include <frogger.hpp>
#include <model_variants.hpp>
#include <serialize_eigen.hpp>

template<typename T, int rank>
Eigen::TensorMap <Eigen::Tensor<T, rank>>
tensor_to_tensor_map(Eigen::Tensor <T, rank> &d) {
  return Eigen::TensorMap < Eigen::Tensor < T, rank
      >> (d.data(), d.dimensions());
}

template<typename ModelVariant>
void save_output(leapfrog::StateSaver<ModelVariant, double> &state_saver,
                 std::string &output_path) {
  auto state = state_saver.get_full_state();

  std::filesystem::path out_path(output_path);
  std::filesystem::path p_total_pop_path = out_path / "p_total_pop";
  serialize::serialize_tensor<double, 3>(state.base.p_total_pop,
                                         p_total_pop_path);

  std::filesystem::path births_path = out_path / "births";
  serialize::serialize_tensor<double, 1>(state.base.births, births_path);

  std::filesystem::path p_total_pop_natural_deaths_path =
      out_path / "p_total_pop_natural_deaths";
  serialize::serialize_tensor<double, 3>(state.base.p_total_pop_natural_deaths,
                                         p_total_pop_natural_deaths_path);

  std::filesystem::path p_hiv_pop_path = out_path / "p_hiv_pop";
  serialize::serialize_tensor<double, 3>(state.base.p_hiv_pop, p_hiv_pop_path);

  std::filesystem::path p_hiv_pop_natural_deaths_path =
      out_path / "p_hiv_pop_natural_deaths";
  serialize::serialize_tensor<double, 3>(state.base.p_hiv_pop_natural_deaths,
                                         p_hiv_pop_natural_deaths_path);

  std::filesystem::path h_hiv_adult_path = out_path / "h_hiv_adult";
  serialize::serialize_tensor<double, 4>(state.base.h_hiv_adult,
                                         h_hiv_adult_path);

  std::filesystem::path h_art_adult_path = out_path / "h_art_adult";
  serialize::serialize_tensor<double, 5>(state.base.h_art_adult,
                                         h_art_adult_path);

  std::filesystem::path h_hiv_deaths_no_art_path =
      out_path / "h_hiv_deaths_no_art";
  serialize::serialize_tensor<double, 4>(state.base.h_hiv_deaths_no_art,
                                         h_hiv_deaths_no_art_path);

  std::filesystem::path p_infections_path = out_path / "p_infections";
  serialize::serialize_tensor<double, 3>(state.base.p_infections,
                                         p_infections_path);

  std::filesystem::path h_hiv_deaths_art_path = out_path / "h_hiv_deaths_art";
  serialize::serialize_tensor<double, 5>(state.base.h_hiv_deaths_art,
                                         h_hiv_deaths_art_path);

  std::filesystem::path h_art_initiation_path = out_path / "h_art_initiation";
  serialize::serialize_tensor<double, 4>(state.base.h_art_initiation,
                                         h_art_initiation_path);

  std::filesystem::path p_hiv_deaths_path = out_path / "p_hiv_deaths";
  serialize::serialize_tensor<double, 3>(state.base.p_hiv_deaths,
                                         p_hiv_deaths_path);
}

int main(int argc, char *argv[]) {
  if (argc < 4) {
    std::cout <<
              "Usage: simulate_model <sim_years> <hts_per_year> <intput_dir> <output_dir>"
              <<
              std::endl;
    return 1;
  }

  int sim_years = atoi(argv[1]);
  int hts_per_year = atoi(argv[2]);
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
      std::cout << "Created output directory '" << output_abs << "'"
                << std::endl;
    } else {
      std::cout << "Failed to create output directory '" << output_abs << "'"
                << std::endl;
    }
  } else {
    std::cout << "Writing to existing output directory " << output_abs << "'"
              << std::endl;
  }

  if (sim_years > 60) {
    std::cout << "Running to max no of sim years: 60\n" << std::endl;
    sim_years = 60;
  }
  if (hts_per_year > 10) {
    std::cout << "Running max no of HIV steps per years: 10" << std::endl;
    hts_per_year = 10;
  }

  // Set working dir to read files
  std::filesystem::path old_dir = std::filesystem::current_path();
  std::filesystem::current_path(input_abs);

  // Only fine-grained ages at first
  const auto ss = leapfrog::StateSpace<leapfrog::BaseModelFullAgeStratification>().base;

  const leapfrog::Options<double> options = {
      hts_per_year,        // HIV steps per year
      30,                  // Time ART start
      ss.hAG,              // Age groups HIV 15+
      // Projection period, 0 for calendar year, 1 for midyear
      leapfrog::internal::PROJPERIOD_CALENDAR
  };

  leapfrog::Tensor1<int> v = serialize::deserialize_tensor<int, 1>(
      std::string("idx_hm_elig"));
  for (int i = 0; i <= sim_years; ++i) {
    // 0-based indexing in C++ vs 1-based indexing in R
    v(i) = v[i] - 1;
  }
  const leapfrog::TensorMap1<int> idx_hm_elig = tensor_to_tensor_map<int, 1>(v);

  leapfrog::Tensor1<double> h(ss.hTS - 1);
  for (int i = 0; i < ss.hTS - 1; ++i) {
    h(i) = 0.5;
  }
  const leapfrog::TensorMap1<double> h_art_stage_dur = tensor_to_tensor_map<double, 1>(
      h);

  leapfrog::Tensor2<double> base_pop_data = serialize::deserialize_tensor<double, 2>(
      std::string("basepop"));
  const leapfrog::TensorMap2<double> base_pop = tensor_to_tensor_map<double, 2>(
      base_pop_data);
  leapfrog::Tensor3<double> survival_probability_data = serialize::deserialize_tensor<double, 3>(
      std::string("survival_probability"));
  const leapfrog::TensorMap3<double> survival_probability = tensor_to_tensor_map<double, 3>(
      survival_probability_data);
  leapfrog::Tensor3<double> net_migration_data = serialize::deserialize_tensor<double, 3>(
      std::string("net_migration"));
  const leapfrog::TensorMap3<double> net_migration = tensor_to_tensor_map<double, 3>(
      net_migration_data);
  leapfrog::Tensor2<double> age_specific_fertility_rate_data = serialize::deserialize_tensor<double, 2>(
      std::string("age_specific_fertility_rate"));
  const leapfrog::TensorMap2<double> age_specific_fertility_rate = tensor_to_tensor_map<double, 2>(
      age_specific_fertility_rate_data);
  leapfrog::Tensor2<double> births_sex_prop_data = serialize::deserialize_tensor<double, 2>(
      std::string("births_sex_prop"));
  const leapfrog::TensorMap2<double> births_sex_prop = tensor_to_tensor_map<double, 2>(
      births_sex_prop_data);
  leapfrog::Tensor1<double> incidence_rate_data = serialize::deserialize_tensor<double, 1>(
      std::string("incidence_rate"));
  const leapfrog::TensorMap1<double> incidence_rate = tensor_to_tensor_map<double, 1>(
      incidence_rate_data);
  leapfrog::Tensor3<double> incidence_age_rate_ratio_data = serialize::deserialize_tensor<double, 3>(
      std::string("incidence_age_rate_ratio"));
  const leapfrog::TensorMap3<double> incidence_age_rate_ratio = tensor_to_tensor_map<double, 3>(
      incidence_age_rate_ratio_data);
  leapfrog::Tensor1<double> incidence_sex_rate_ratio_data = serialize::deserialize_tensor<double, 1>(
      std::string("incidence_sex_rate_ratio"));
  const leapfrog::TensorMap1<double> incidence_sex_rate_ratio = tensor_to_tensor_map<double, 1>(
      incidence_sex_rate_ratio_data);
  leapfrog::Tensor3<double> cd4_mortality_data = serialize::deserialize_tensor<double, 3>(
      std::string("cd4_mortality_full"));
  const leapfrog::TensorMap3<double> cd4_mortality = tensor_to_tensor_map<double, 3>(
      cd4_mortality_data);
  leapfrog::Tensor3<double> cd4_progression_data = serialize::deserialize_tensor<double, 3>(
      std::string("cd4_progression_full"));
  const leapfrog::TensorMap3<double> cd4_progression = tensor_to_tensor_map<double, 3>(
      cd4_progression_data);
  leapfrog::Tensor3<double> cd4_initial_distribution_data = serialize::deserialize_tensor<double, 3>(
      std::string("cd4_initial_distribution_full"));
  const leapfrog::TensorMap3<double> cd4_initial_distribution = tensor_to_tensor_map<double, 3>(
      cd4_initial_distribution_data);
  leapfrog::Tensor4<double> art_mortality_rate_data = serialize::deserialize_tensor<double, 4>(
      std::string("art_mortality_rate_full"));
  const leapfrog::TensorMap4<double> art_mortality_rate = tensor_to_tensor_map<double, 4>(
      art_mortality_rate_data);
  leapfrog::Tensor2<double> art_mortality_time_rate_ratio_data = serialize::deserialize_tensor<double, 2>(
      std::string("art_mortality_time_rate_ratio"));
  const leapfrog::TensorMap2<double> art_mortality_time_rate_ratio = tensor_to_tensor_map<double, 2>(
      art_mortality_time_rate_ratio_data);
  leapfrog::Tensor1<double> art_dropout_data = serialize::deserialize_tensor<double, 1>(
      std::string("art_dropout"));
  const leapfrog::TensorMap1<double> art_dropout = tensor_to_tensor_map<double, 1>(
      art_dropout_data);
  Eigen::Tensor<double, 2> adults_on_art_data = serialize::deserialize_tensor<double, 2>(
      std::string("adults_on_art"));
  const leapfrog::TensorMap2<double> adults_on_art = tensor_to_tensor_map<double, 2>(
      adults_on_art_data);
  leapfrog::Tensor2<int> adults_on_art_is_percent_data = serialize::deserialize_tensor<int, 2>(
      std::string("adults_on_art_is_percent"));
  const leapfrog::TensorMap2<int> adults_on_art_is_percent = tensor_to_tensor_map<int, 2>(
      adults_on_art_is_percent_data);

  const leapfrog::Demography<double> demography = {
      base_pop,
      survival_probability,
      net_migration,
      age_specific_fertility_rate,
      births_sex_prop
  };

  const leapfrog::Incidence<double> incidence = {
      incidence_rate,
      incidence_age_rate_ratio,
      incidence_sex_rate_ratio
  };

  const leapfrog::NaturalHistory<double> natural_history = {
      cd4_mortality,
      cd4_progression,
      cd4_initial_distribution,
      1   // Scale CD4 mortality
  };

  const leapfrog::Art<double> art = {
      idx_hm_elig,
      art_mortality_rate,
      art_mortality_time_rate_ratio,
      art_dropout,
      adults_on_art,
      adults_on_art_is_percent,
      h_art_stage_dur,
      0.2  // initiation_mortality_weight
  };

  const leapfrog::Parameters<leapfrog::BaseModelFullAgeStratification, double> params = {
      options,
      demography,
      incidence,
      natural_history,
      art};

  leapfrog::internal::IntermediateData<leapfrog::BaseModelFullAgeStratification, double> intermediate(
      options.hAG_15plus);

  std::vector<int> save_steps(61);
  std::iota(save_steps.begin(), save_steps.end(), 0);
  leapfrog::StateSaver<leapfrog::BaseModelFullAgeStratification, double> state_output(
      sim_years, save_steps);

  const char *n_runs_char = std::getenv("N_RUNS");
  size_t n_runs = 1;
  if (n_runs_char != nullptr) {
    // If we're profiling we want to get accurate info about where time is spent during the
    // main model fit. This runs so quickly though that just going through once won't sample enough
    // times for us to see. And it will sample from the tensor file serialization/deserialization more.
    // So we run the actual model fit multiple times when profiling so the sampler can actually pick
    // up the slow bits.
    n_runs = atoi(n_runs_char);
    std::cout << "Running model fit " << n_runs << " times" << std::endl;
  }

  for (size_t i = 0; i < n_runs; ++i) {
    auto state_current = leapfrog::State<leapfrog::BaseModelFullAgeStratification, double>(
        params);
    auto state_next = state_current;
    leapfrog::set_initial_state<leapfrog::BaseModelFullAgeStratification, double>(state_current, params);

    // Save initial state
    state_output.save_state(state_current, 0);

    // Each time step is mid-point of the year
    for (int step = 1; step <= sim_years; ++step) {
      leapfrog::run_general_pop_demographic_projection<leapfrog::BaseModelFullAgeStratification, double>(
          step, params,
          state_current,
          state_next,
          intermediate);
      leapfrog::run_hiv_pop_demographic_projection<leapfrog::BaseModelFullAgeStratification, double>(
          step, params,
          state_current,
          state_next,
          intermediate);
      leapfrog::run_hiv_model_simulation<leapfrog::BaseModelFullAgeStratification, double>(
          step, params, state_current,
          state_next, intermediate);
      state_output.save_state(state_next, step);
      std::swap(state_current, state_next);
      intermediate.reset();
      state_next.reset();
    }
  }
  std::cout << "Fit complete" << std::endl;

  save_output(state_output, output_abs);

  return 0;
}
