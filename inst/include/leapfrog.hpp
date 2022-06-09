#pragma once

#include <consts.hpp>
#include <typedefs.hpp>
#include <unsupported/Eigen/CXX11/Tensor>

using Eigen::Sizes;
using Eigen::Tensor;
using Eigen::TensorBase;
using Eigen::TensorFixedSize;
using Eigen::TensorMap;

inline int add(int a, int b) {
  return a + b;
}

//' The output memory is passed as an argument rather than constructed
//' by the function.
//'
//' Template parameters
//' @param Type data type variables involved in calculations (typically 'double'
//'   but required to be a templated parameter for autodiff integration).
//' @param NUM_GENDERS number of genders (= 2)
//' @param AGE_GROUPS_POP number of population age groups (e.g. 81 for ages 0 to
//'    80+)
//' @param FERTILITY_FIRST_AGE_GROUP first index eligible for fertility
//' @param AGE_GROUPS_FERT number of ages eligible for fertility
//'
//' @details
//' State space dimensions are specified as template parameters so that
//' these are specified at compile time, allowing stack allocation of
//' working arrays.
//'
//' Notes:
//' * Consider adding half the migrants at the start and half at the end?
//' * Unsure whether to have this function accept raw pointers or TensorMap
//'   - If passing Tensor objects, it would be ideal to pass TensorBase, but
//'     this is not possible because dimensions are not known.

template <typename Type,
          int NUM_GENDERS,
          int AGE_GROUPS_POP,
          int FERTILITY_FIRST_AGE_GROUP,
          int AGE_GROUPS_FERT>
void leapfrog_sim(const Type* p_base_pop,
                  const Type* p_survival,
                  const Type* p_net_migration,
                  const Type* p_asfr,
                  const Type* p_births_sex_prop,
                  //
                  // settings
                  const std::string model_type,
                  const int sim_years,
                  //
                  // outputs
                  Type* p_total_population,
                  Type* p_births,
                  Type* p_natural_deaths) {
  Model model;
  if (model_type == "base") {
    model = BaseModel();
  } else {
    throw Exception  // TODO: make this a better exception
  }
  model.run_model()
}
}
// // inputs

// demography
const TensorMapX2cT base_pop(p_base_pop, AGE_GROUPS_POP, NUM_GENDERS);
const TensorMapX3cT survival(p_survival,
                             AGE_GROUPS_POP + 1,
                             NUM_GENDERS,
                             sim_years);
const TensorMapX3cT net_migration(p_net_migration,
                                  AGE_GROUPS_POP,
                                  NUM_GENDERS,
                                  sim_years);
const TensorMapX2cT asfr(p_asfr, AGE_GROUPS_FERT, sim_years);
const TensorMapX2cT births_sex_prop(p_births_sex_prop, NUM_GENDERS, sim_years);

// outputs
TensorMapX3T total_population(p_total_population,
                              AGE_GROUPS_POP,
                              NUM_GENDERS,
                              sim_years);
TensorMapX1T births(p_births, sim_years);
TensorMapX3T natural_deaths(p_natural_deaths,
                            AGE_GROUPS_POP,
                            NUM_GENDERS,
                            sim_years);

// initialise population

for (int g = 0; g < NUM_GENDERS; g++) {
  for (int a = 0; a < AGE_GROUPS_POP; a++) {
    total_population(a, g, 0) = base_pop(a, g);
  }
}
}

void run_demographic_projection() {
  TensorFixedSize<Type, Sizes<AGE_GROUPS_POP, NUM_GENDERS>> migrate_ag;

  // ageing and non-HIV mortality
  for (int g = 0; g < NUM_GENDERS; g++) {
    for (int a = 1; a < AGE_GROUPS_POP; a++) {
      natural_deaths(a, g, t) =
          total_population(a - 1, g, t - 1) * (1.0 - survival(a, g, t));
      total_population(a, g, t) =
          total_population(a - 1, g, t - 1) - natural_deaths(a, g, t);
    }

    // open age group
    Type natural_deaths_open_age =
        total_population(AGE_GROUPS_POP - 1, g, t - 1) *
        (1.0 - survival(AGE_GROUPS_POP, g, t));
    natural_deaths(AGE_GROUPS_POP - 1, g, t) += natural_deaths_open_age;
    total_population(AGE_GROUPS_POP - 1, g, t) +=
        total_population(AGE_GROUPS_POP - 1, g, t - 1) -
        natural_deaths_open_age;

    // net migration
    for (int a = 1; a < AGE_GROUPS_POP - 1; a++) {
      // Number of net migrants adjusted for survivorship to end of period (qx
      // / 2)
      migrate_ag(a, g) = net_migration(a, g, t) * (1.0 + survival(a, g, t)) *
                         0.5 / total_population(a, g, t);
      total_population(a, g, t) *= 1.0 + migrate_ag(a, g);
    }

    // For open age group, net_migrationant survivor adjustment based on
    // weighted survival for age 79 and age 80+.
    // * Numerator: total_population(a, g, t-1) * (1.0 + survival(a+1, g, t))
    // + total_population(a-1, g, t-1) * (1.0 + survival(a, g, t))
    // * Denominator: total_population(a, g, t-1) + total_population(a-1, g,
    // t-1) Re-expressed current population and deaths to open age group
    // (already calculated):
    int a = AGE_GROUPS_POP - 1;
    Type survival_netmig =
        (total_population(a, g, t) +
         0.5 * natural_deaths(AGE_GROUPS_POP - 1, g, t)) /
        (total_population(a, g, t) + natural_deaths(AGE_GROUPS_POP - 1, g, t));
    migrate_ag(a, g) =
        survival_netmig * net_migration(a, g, t) / total_population(a, g, t);
    total_population(a, g, t) *= 1.0 + migrate_ag(a, g);
  }

  // fertility

  births(t) = 0.0;
  for (int af = 0; af < AGE_GROUPS_FERT; af++) {
    births(t) +=
        (total_population(FERTILITY_FIRST_AGE_GROUP + af, FEMALE, t - 1) +
         total_population(FERTILITY_FIRST_AGE_GROUP + af, FEMALE, t)) *
        0.5 * asfr(af, t);
  }

  // add births
  for (int g = 0; g < NUM_GENDERS; g++) {
    Type births_sex = births(t) * births_sex_prop(g, t);
    natural_deaths(0, g, t) = births_sex * (1.0 - survival(0, g, t));
    total_population(0, g, t) = births_sex * survival(0, g, t);

    // Assume 2/3 survival rate since mortality in first six months higher
    // than second 6 months (Spectrum manual, section 6.2.7.4)
    Type migrate_a0 = net_migration(0, g, t) * (1.0 + 2.0 * survival(0, g, t)) /
                      3.0 / total_population(0, g, t);
    total_population(0, g, t) *= 1.0 + migrate_a0;
  }
}