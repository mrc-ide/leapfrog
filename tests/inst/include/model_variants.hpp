#pragma once

namespace leapfrog {

enum HivAgeStratification {
    full,
    coarse
};

struct BaseModelFullAgeStratification {
  static constexpr HivAgeStratification stratification = full;
  static constexpr bool run_child_model = false;
};

struct BaseModelCoarseAgeStratification {
  static constexpr HivAgeStratification stratification = coarse;
  static constexpr bool run_child_model = false;
};

struct ChildModel {
  static constexpr HivAgeStratification stratification = full;
  static constexpr bool run_child_model = true;
};

}
