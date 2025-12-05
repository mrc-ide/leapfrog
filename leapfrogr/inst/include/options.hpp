#pragma once

#include <algorithm>
#include <stdexcept>
#include <string_view>
#include <vector>
#include <cmath>
#include <string>

#include "generated/state_space.hpp"

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

int get_proj_period_enum(std::string_view projection_period) {
  if (projection_period == "midyear") {
    return BaseSS::PROJPERIOD_MIDYEAR;
  } else if (projection_period == "calendar") {
    return BaseSS::PROJPERIOD_CALENDAR;
  } else {
    throw std::invalid_argument(
      "Invalid projection period: '" + std::string(projection_period) +
      "'. Allowed values are: 'midyear' or 'calendar'.");
  }
}

}

template<typename real_type>
struct Options {
  int hts_per_year;
  double dt;
  const int ts_art_start;
  const int proj_period_int;
  const int proj_start_year;
  const int proj_end_year;
  const int proj_steps;

  Options(
    int hts_per_year,
    int ts_art_start,
    int proj_period_int,
    int proj_start_year,
    int proj_end_year
  ):
    hts_per_year(hts_per_year),
    dt(1.0 / hts_per_year),
    ts_art_start(ts_art_start),
    proj_period_int(proj_period_int),
    proj_start_year(proj_start_year),
    proj_end_year(proj_end_year),
    proj_steps(proj_end_year - proj_start_year + 1) {}
};

template<typename real_type>
const Options<real_type> get_opts(
  const int hiv_steps,
  const int t_ART_start,
  const std::string_view projection_period,
  const int proj_start_year,
  const std::vector<int>& output_years
) {
  internal::validate_output_years(output_years, proj_start_year);
  const int proj_period = internal::get_proj_period_enum(projection_period);
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
