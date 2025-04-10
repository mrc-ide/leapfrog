#pragma once

#include "generated/state_space_mixer.hpp"

namespace leapfrog {

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

  Options(
    int hts_per_year,
    int ts_art_start,
    int proj_period_int
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
    proj_period_int(proj_period_int) {}
};

template<typename real_type>
const Options<real_type> get_opts(
  const int hiv_steps,
  const int t_ART_start,
  const bool is_midyear_projection
) {
  const int proj_period = is_midyear_projection
    ? internal::BaseSS::PROJPERIOD_MIDYEAR
    : internal::BaseSS::PROJPERIOD_CALENDAR;
  const Options<real_type> opts = {
    hiv_steps,
    t_ART_start,
    proj_period
  };
  return opts;
};

}
