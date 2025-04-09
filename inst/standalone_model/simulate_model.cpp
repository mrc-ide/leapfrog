#include <cstdlib>
#include <filesystem>
#include <vector>

#include "frogger.hpp"
#include "options.hpp"
#include "generated/model_variants.hpp"
#include "generated/cpp_interface/cpp_adapter_mixer.hpp"
#include "serialize_eigen.hpp"

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
    std::cout << "Input dir '" << input_dir << "' does not exist." << std::endl;
    return 1;
  }

  std::filesystem::path output_abs = std::filesystem::absolute(output_dir);
  if (!std::filesystem::exists(output_abs)) {
    if (std::filesystem::create_directory(output_abs)) {
      std::cout << "Created output directory '" << std::string{output_abs} << "'"
                << std::endl;
    } else {
      std::cout << "Failed to create output directory '" << std::string{output_abs} << "'"
                << std::endl;
    }
  } else {
    std::cout << "Writing to existing output directory " << std::string{output_abs} << "'"
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
  std::iota(save_steps.begin(), save_steps.end(), 0);

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

  const auto opts = leapfrog::get_opts<double>(hts_per_year, 30, true);
  const auto pars = leapfrog::AdapterCpp<double, leapfrog::HivFullAgeStratification>::get_pars(input_dir, opts, sim_years);

  for (size_t i = 0; i < n_runs; ++i) {
    auto state = leapfrog::Leapfrog<double, leapfrog::HivFullAgeStratification>::run_model(sim_years, save_steps, pars, opts);
  }
  std::cout << "Fit complete" << std::endl;

  auto state = leapfrog::Leapfrog<double, leapfrog::HivFullAgeStratification>::run_model(sim_years, save_steps, pars, opts);
  leapfrog::AdapterCpp<double, leapfrog::HivFullAgeStratification>::build_output(0, state, output_abs);

  return 0;
}
