#pragma once

namespace leapfrog {

struct DemographicProjection {
  static constexpr bool run_demographic_projection = true;
  static constexpr bool run_hiv_simulation = false;
  static constexpr bool use_coarse_stratification = false;
  static constexpr bool run_child_model = false;
};
struct HivFullAgeStratification {
  static constexpr bool run_demographic_projection = true;
  static constexpr bool run_hiv_simulation = true;
  static constexpr bool use_coarse_stratification = false;
  static constexpr bool run_child_model = false;
};
struct HivCoarseAgeStratification {
  static constexpr bool run_demographic_projection = true;
  static constexpr bool run_hiv_simulation = true;
  static constexpr bool use_coarse_stratification = true;
  static constexpr bool run_child_model = false;
};
struct ChildModel {
  static constexpr bool run_demographic_projection = true;
  static constexpr bool run_hiv_simulation = true;
  static constexpr bool use_coarse_stratification = false;
  static constexpr bool run_child_model = true;
};

}