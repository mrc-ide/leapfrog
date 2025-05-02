#pragma once

namespace leapfrog {
namespace internal {

#pragma pack(push, 8)
struct DpParams {
  double* base_pop;
  int base_pop_length;
  double* survival_probability;
  int survival_probability_length;
  double* net_migration;
  int net_migration_length;
  double* age_specific_fertility_rate;
  int age_specific_fertility_rate_length;
  double* births_sex_prop;
  int births_sex_prop_length;
};

struct DpOut {
  double* p_total_pop;
  int p_total_pop_length;
  double* p_total_pop_natural_deaths;
  int p_total_pop_natural_deaths_length;
  double* births;
  int births_length;
};
#pragma pack(pop)

}
}
