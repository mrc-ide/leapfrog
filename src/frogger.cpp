#include <Rcpp.h>

#include "frogger.hpp"
#include "types.hpp"
#include "state_space.hpp"
#include "model_setup.hpp"
#include "model_output.hpp"

int transform_simulation_years(const Rcpp::List demp, SEXP r_sim_years) {
  Rcpp::NumericVector Sx = demp["Sx"];
  Rcpp::Dimension d = Sx.attr("dim");
  // Simulation initialises state from first years input data (index 0)
  // then runs for each year simulating this years (i) data using previous years
  // state (i - 1) and this years input data (i). So -1 off index for max years
  // to simulate as index 0 used for initial state
  const int max_sim_years = d[2] - 1;
  if (r_sim_years == R_NilValue) {
    return max_sim_years;
  }
  auto sim_years = INTEGER(r_sim_years)[0];
  if (sim_years > max_sim_years) {
    Rcpp::stop("No of years > max years of " + std::to_string(max_sim_years));
  }
  return sim_years;
}

int transform_hiv_steps_per_year(SEXP r_hiv_steps_per_year) {
  int hiv_steps_per_year;
  if (r_hiv_steps_per_year == R_NilValue) {
    hiv_steps_per_year = 10;
  } else {
    hiv_steps_per_year = INTEGER(r_hiv_steps_per_year)[0];
  }
  return hiv_steps_per_year;
}

template <leapfrog::HivAgeStratification S>
leapfrog::StateSpace<S> get_state_space(const std::string hiv_age_stratification) {
  leapfrog::StateSpace<S> ss;
  if (hiv_age_stratification == "full") {
    ss = leapfrog::StateSpace<leapfrog::HivAgeStratification::full>();
  } else {
    // We've already validated that hiv_age_stratification is one of "full" or "coarse"
    ss = leapfrog::StateSpace<leapfrog::HivAgeStratification::coarse>();
  }
  return ss;
}

std::vector<int> transform_output_steps(Rcpp::NumericVector output_steps) {
  return Rcpp::as<std::vector<int>>(output_steps);
}

template <leapfrog::HivAgeStratification S>
Rcpp::List fit_model(const leapfrog::StateSpace<S> ss,
                     const Rcpp::List data,
                     const int proj_years,
                     const int hiv_steps,
                     const std::vector<int> save_steps) {
  const leapfrog::Options<double> opts = {
      hiv_steps,
      Rcpp::as<int>(data["t_ART_start"]) - 1,
      ss.age_groups_hiv,
      Rcpp::as<int>(data["scale_cd4_mort"]),
      Rcpp::as<double>(data["art_alloc_mxweight"])
  };

  const leapfrog::Parameters<double> params = setup_model_params<double, S>(data, opts, proj_years);

  auto state = leapfrog::run_model<double, S>(proj_years, save_steps, params);

  auto ret = build_r_output<double, S>(state, save_steps);

  return ret;
}

// [[Rcpp::export]]
Rcpp::List run_base_model(const Rcpp::List data,
                          SEXP sim_years,
                          SEXP hiv_steps_per_year,
                          Rcpp::NumericVector output_steps,
                          std::string hiv_age_stratification) {
  const int proj_years = transform_simulation_years(data, sim_years);
  const std::vector<int> save_steps = transform_output_steps(output_steps);
  const int hiv_steps = transform_hiv_steps_per_year(hiv_steps_per_year);


  Rcpp::List ret;
  if (hiv_age_stratification == "full") {
    constexpr auto ss = leapfrog::StateSpace<leapfrog::HivAgeStratification::full>();
    ret = fit_model<leapfrog::HivAgeStratification::full>(ss, data,
                                                          proj_years, hiv_steps, save_steps);
  } else {
    // We've already validated that hiv_age_stratification is one of "full" or "coarse"
    constexpr auto ss = leapfrog::StateSpace<leapfrog::HivAgeStratification::coarse>();
    ret = fit_model<leapfrog::HivAgeStratification::coarse>(ss, data,
                                                            proj_years, hiv_steps, save_steps);
  }

  return ret;
}
