#pragma once

#include <algorithm>

#include "generated/state_space_mixer.hpp"

namespace leapfrog {

namespace internal {

void validate_output_years(const std::vector<int> output_years,
                                 const int proj_start_year) {
  const auto first_year = std::min_element(std::begin(output_years),
                                           std::end(output_years));
  if (*first_year < proj_start_year) {
    throw std::runtime_error("Trying to output for year: '" +
        std::to_string(*first_year) +
        "' which is before the projection start year: '" +
        std::to_string(proj_start_year) +
        "'.");
  }
}

int get_proj_end_year(const std::vector<int>& output_years) {
  const auto last_year = std::max_element(std::begin(output_years),
                                          std::end(output_years));
  return *last_year;
}

}

template<typename real_type>
struct Options {
  int hts_per_year;
  double dt;
  const int p_idx_fertility_first;
  const int p_fertility_age_groups;
  const int p_idx_hiv_first_adult;
  const int adult_incidence_first_age_group;
  const int pAG_INCIDPOP;
  const int ts_art_start;
  const int hIDX_15PLUS;
  const int proj_period_int;
  const int proj_start_year;
  const int proj_end_year;
  const int proj_time_steps;

  Options(
    int hts_per_year,
    int ts_art_start,
    int proj_period_int,
    int proj_start_year,
    int proj_end_year
  ):
    hts_per_year(hts_per_year),
    dt(1.0 / hts_per_year),
    p_idx_fertility_first(15),
    p_fertility_age_groups(35),
    p_idx_hiv_first_adult(15),
    adult_incidence_first_age_group(15),
    pAG_INCIDPOP(35),
    ts_art_start(ts_art_start),
    hIDX_15PLUS(0),
    proj_period_int(proj_period_int),
    proj_start_year(proj_start_year),
    proj_end_year(proj_end_year),
    proj_time_steps(proj_end_year - proj_start_year + 1) {}
};

template<typename real_type>
const Options<real_type> get_opts(
  const int hiv_steps,
  const int t_ART_start,
  const bool is_midyear_projection,
  const int proj_start_year,
  const std::vector<int>& output_years
) {
  internal::validate_output_years(output_years, proj_start_year);
  const int proj_period = is_midyear_projection
    ? internal::BaseSS::PROJPERIOD_MIDYEAR
    : internal::BaseSS::PROJPERIOD_CALENDAR;
  const int proj_end_year = internal::get_proj_end_year(output_years);
  const Options<real_type> opts = {
    hiv_steps,
    t_ART_start,
    proj_period,
    proj_start_year,
    proj_end_year
  };
  return opts;
};

} // namespace leapfrog
