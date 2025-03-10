#pragma once

namespace leapfrog {

enum HivAgeStratification {
    full,
    coarse,
    none
};

struct DemographicProjection {
  static constexpr bool run_hiv_simulation = false;
  static constexpr HivAgeStratification stratification = none;
  static constexpr bool run_child_model = false;
};

struct HivFullAgeStratification {
  static constexpr bool run_hiv_simulation = true;
  static constexpr HivAgeStratification stratification = full;
  static constexpr bool run_child_model = false;
};

struct HivCoarseAgeStratification {
  static constexpr bool run_hiv_simulation = true;
  static constexpr HivAgeStratification stratification = coarse;
  static constexpr bool run_child_model = false;
};

struct ChildModel {
  static constexpr bool run_hiv_simulation = true;
  static constexpr HivAgeStratification stratification = full;
  static constexpr bool run_child_model = true;
};

}
