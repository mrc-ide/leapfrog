#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <vector>

#include "frogger.hpp"
#include "serialize_eigen.hpp"
#include "generated/config_mixer.hpp"

template<typename T, int rank>
Eigen::TensorMap <Eigen::Tensor<T, rank>>
tensor_to_tensor_map(Eigen::Tensor <T, rank> &d) {
  return Eigen::TensorMap < Eigen::Tensor < T, rank
      >> (d.data(), d.dimensions());
}

template<typename real_type, typename ModelVariant>
void save_output(const typename leapfrog::ConfigMixed<real_type, ModelVariant>::OutputState &state,
                 std::string &output_path) {
  std::filesystem::path out_path(output_path);
  std::filesystem::path p_total_pop_path = out_path / "p_total_pop";
  serialize::serialize_tensor<double, 3>(state.p_total_pop,
                                         p_total_pop_path);

  std::filesystem::path births_path = out_path / "births";
  serialize::serialize_tensor<double, 1>(state.births, births_path);

  std::filesystem::path p_total_pop_natural_deaths_path =
      out_path / "p_total_pop_natural_deaths";
  serialize::serialize_tensor<double, 3>(state.p_total_pop_natural_deaths,
                                         p_total_pop_natural_deaths_path);

  std::filesystem::path p_hiv_pop_path = out_path / "p_hiv_pop";
  serialize::serialize_tensor<double, 3>(state.p_hiv_pop, p_hiv_pop_path);

  std::filesystem::path p_hiv_pop_natural_deaths_path =
      out_path / "p_hiv_pop_natural_deaths";
  serialize::serialize_tensor<double, 3>(state.p_hiv_pop_natural_deaths,
                                         p_hiv_pop_natural_deaths_path);

  std::filesystem::path h_hiv_adult_path = out_path / "h_hiv_adult";
  serialize::serialize_tensor<double, 4>(state.h_hiv_adult,
                                         h_hiv_adult_path);

  std::filesystem::path h_art_adult_path = out_path / "h_art_adult";
  serialize::serialize_tensor<double, 5>(state.h_art_adult,
                                         h_art_adult_path);

  std::filesystem::path h_hiv_deaths_no_art_path =
      out_path / "h_hiv_deaths_no_art";
  serialize::serialize_tensor<double, 4>(state.h_hiv_deaths_no_art,
                                         h_hiv_deaths_no_art_path);

  std::filesystem::path p_infections_path = out_path / "p_infections";
  serialize::serialize_tensor<double, 3>(state.p_infections,
                                         p_infections_path);

  std::filesystem::path h_hiv_deaths_art_path = out_path / "h_hiv_deaths_art";
  serialize::serialize_tensor<double, 5>(state.h_hiv_deaths_art,
                                         h_hiv_deaths_art_path);

  std::filesystem::path h_art_initiation_path = out_path / "h_art_initiation";
  serialize::serialize_tensor<double, 4>(state.h_art_initiation,
                                         h_art_initiation_path);

  std::filesystem::path p_hiv_deaths_path = out_path / "p_hiv_deaths";
  serialize::serialize_tensor<double, 3>(state.p_hiv_deaths,
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

  std::filesystem::path input_abs = std::filesystem::absolute(input_dir);
  if (!std::filesystem::exists(input_abs)) {
    std::cout << "Input dir '" << input_dir << "' does not exist.\n";
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

  if (sim_years > 61) {
    std::cout << "Running to max no of sim years: 61\n" << std::endl;
    sim_years = 61;
  }
  if (hts_per_year > 10) {
    std::cout << "Running max no of HIV steps per years: 10" << std::endl;
    hts_per_year = 10;
  }

  std::vector<int> save_steps(sim_years);
  const auto state = leapfrog::Leapfrog<double, leapfrog::HivFullAgeStratification>::simulate_model(input_abs, sim_years, hts_per_year, save_steps, true, 30);

  save_output<double, leapfrog::HivFullAgeStratification>(state, output_abs);

  return 0;
}
