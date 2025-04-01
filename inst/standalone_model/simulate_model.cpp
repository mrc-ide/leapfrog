#include <filesystem>
#include <fstream>
#include <vector>

#include "frogger.hpp"
#include "generated/config_mixer.hpp"

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
  std::iota(save_steps.begin(), save_steps.end(), 0);

  leapfrog::Leapfrog<double, leapfrog::HivFullAgeStratification>::simulate_model(input_abs, sim_years, hts_per_year, save_steps, true, 30, output_abs);

  std::cout << "Fit complete" << std::endl;

  return 0;
}
