#pragma once

namespace leapfrog {
namespace internal {

#pragma pack(push, 8)
struct DpParams {
    double* base_pop;
    double* survival_probability;
    double* net_migration;
    double* age_specific_fertility_rate;
    double* births_sex_prop;
};

struct DpOut {
    double* p_total_pop;
    double* p_total_pop_natural_deaths;
    double* births;
};
#pragma pack(pop)

}
}
