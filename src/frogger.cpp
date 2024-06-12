#include <Rcpp.h>

#include "frogger.hpp"
#include "types.hpp"
#include "model_variants.hpp"
#include "state_space.hpp"
#include "model_input.hpp"
#include "model_output.hpp"
#include "r_utils.hpp"

int transform_simulation_years(
    const Rcpp::List demp,
    Rcpp::Nullable<Rcpp::NumericVector> r_sim_years) {

  Rcpp::NumericVector Sx = demp["Sx"];
  Rcpp::Dimension d = Sx.attr("dim");
  // Simulation initialises state from first years input data (index 0)
  // then runs for each year simulating this years (i) data using previous years
  // state (i - 1) and this years input data (i).
  const int max_sim_years = d[2];
  if (r_sim_years.isNull()) {
    return max_sim_years;
  }
  Rcpp::NumericVector years(r_sim_years);
  auto sim_years = years.length();
  if (sim_years > max_sim_years) {
    Rcpp::stop("No of years > max years of " + std::to_string(max_sim_years));
  }
  return sim_years;
}

int transform_hts_per_year(SEXP r_hts_per_year) {
  int hts_per_year;
  if (r_hts_per_year == R_NilValue) {
    hts_per_year = 10;
  } else {
    hts_per_year = INTEGER(r_hts_per_year)[0];
  }
  return hts_per_year;
}

std::vector<int> transform_output_steps(Rcpp::NumericVector output_steps) {
  auto out = Rcpp::as<std::vector<int>>(output_steps);
  convert_0_based(out);
  return out;
}

template<typename ModelVariant>
Rcpp::List simulate_model(const leapfrog::StateSpace<ModelVariant> ss,
                          const Rcpp::List data,
                          const int proj_years,
                          const int hiv_steps,
                          const std::vector<int> save_steps) {
  const leapfrog::Options<double> opts = {
      hiv_steps,
      Rcpp::as<int>(data["t_ART_start"]) - 1,
      ss.base.hAG
  };

  const auto params = setup_model_params<ModelVariant, double>(data, opts, proj_years);

  auto state = leapfrog::run_model<ModelVariant, double>(proj_years, save_steps, params);

  auto ret = build_r_output<ModelVariant, double>(state, save_steps);

  return ret;
}

// [[Rcpp::export]]
Rcpp::List run_base_model(const Rcpp::List data,
                          Rcpp::Nullable<Rcpp::NumericVector> sim_years,
                          SEXP hts_per_year,
                          Rcpp::NumericVector output_steps,
                          std::string model_variant) {
  const int proj_years = transform_simulation_years(data, sim_years);
  const std::vector<int> save_steps = transform_output_steps(output_steps);
  const int hiv_steps = transform_hts_per_year(hts_per_year);

  Rcpp::List ret;
  if (model_variant == "ChildModel") {
    constexpr auto ss = leapfrog::StateSpace<leapfrog::ChildModel>();
    ret = simulate_model<leapfrog::ChildModel>(ss, data, proj_years, hiv_steps, save_steps);
  } else if (model_variant == "BaseModelFullAgeStratification") {
    constexpr auto ss = leapfrog::StateSpace<leapfrog::BaseModelFullAgeStratification>();
    ret = simulate_model<leapfrog::BaseModelFullAgeStratification>(ss, data, proj_years, hiv_steps, save_steps);
  } else if (model_variant == "BaseModelCoarseAgeStratification") {
    constexpr auto ss = leapfrog::StateSpace<leapfrog::BaseModelCoarseAgeStratification>();
    ret = simulate_model<leapfrog::BaseModelCoarseAgeStratification>(ss, data, proj_years, hiv_steps, save_steps);
  } else {
    Rcpp::stop("Invalid model variant " + model_variant + " must be one of " +
               "'BaseModelFullAgeStratification', 'BaseModelCoarseAgeStratification' or 'ChildModel'");
  }

  return ret;
}

// [[Rcpp::export]]
Rcpp::List run_wlhiv_births_r(const Rcpp::List data,
                              const Rcpp::Nullable<Rcpp::NumericVector> sim_years,
                              const SEXP hts_per_year,
                              const Rcpp::NumericVector output_steps,
                              const Rcpp::List data_from_adult_model) {

    constexpr auto ss = leapfrog::StateSpace<leapfrog::ChildModel>();
    const int proj_years = transform_simulation_years(data, sim_years);
    const std::vector<int> save_steps= transform_output_steps(output_steps);
    const int hiv_steps = transform_hts_per_year(hts_per_year);
    
    const leapfrog::Options<double> opts = {
        hiv_steps,
        Rcpp::as<int>(data["t_ART_start"]) - 1,
        ss.base.hAG
    };

    const auto pars = setup_model_params<leapfrog::ChildModel, double>(data, opts, proj_years);

    size_t output_years = save_steps.size();
    auto state = leapfrog::State<leapfrog::ChildModel, double>(pars);

    const leapfrog::TensorMap4<double> h_hiv_adult = parse_data<double>(data_from_adult_model, "h_hiv_adult", ss.base.hDS, ss.base.hAG, ss.base.NS, output_years);
    const leapfrog::TensorMap5<double> h_art_adult = parse_data<double>(data_from_adult_model, "h_art_adult", ss.base.hTS, ss.base.hDS, ss.base.hAG, ss.base.NS, output_years);
    const leapfrog::TensorMap3<double> p_total_pop = parse_data<double>(data_from_adult_model, "p_total_pop", ss.base.pAG, ss.base.NS, output_years);
    const leapfrog::TensorMap1<double> births = parse_data<double>(data_from_adult_model, "births", output_years);

    auto state_next = state;
    leapfrog::set_initial_state<leapfrog::ChildModel, double>(state, pars);
    leapfrog::internal::IntermediateData<leapfrog::ChildModel, double> intermediate(pars.base.options.hAG_15plus);
    intermediate.reset();

    leapfrog::StateSaver<leapfrog::ChildModel, double> state_output(proj_years, save_steps);
    // Save initial state
    state_output.save_state(state, 0);

    // Each time step is mid-point of the year
    for (int step = 1; step < proj_years; ++step) {

      // sub in values from our previous adult run
      for (int hd = 0; hd < ss.base.hDS; ++hd) {
        for (int a = 0; a < ss.base.hAG; ++a) {
          for (int s = 0; s < ss.base.NS; ++s) {
            state_next.base.h_hiv_adult(hd, a, s) = h_hiv_adult(hd, a, s, step);
            state.base.h_hiv_adult(hd, a, s) = h_hiv_adult(hd, a, s, step - 1);
          }
        }
      }
      for (int hu = 0; hu < ss.base.hTS; ++hu) {
        for (int hd = 0; hd < ss.base.hDS; ++hd) {
          for (int a = 0; a < ss.base.hAG; ++a) {
            for (int s = 0; s < ss.base.NS; ++s) {
              state_next.base.h_art_adult(hu, hd, a, s) = h_art_adult(hu, hd, a, s, step);
              state.base.h_art_adult(hu, hd, a, s) = h_art_adult(hu, hd, a, s, step - 1);
            }
          }
        }
      }
      for (int p = 0; p < ss.base.pAG; ++p) {
        for (int s = 0; s < ss.base.NS; ++s) {
          state_next.base.p_total_pop(p, s) = p_total_pop(p, s, step);
          state.base.p_total_pop(p, s) = p_total_pop(p, s, step - 1);
        }
      }
      state_next.base.births = births(step);
      state.base.births = births(step - 1);

      run_wlhiv_births(step, pars, state, state_next, intermediate);

      state_output.save_state(state_next, step);
      std::swap(state, state_next);
      state_next.reset();
      intermediate.reset();
    }

    auto full_state_output = state_output.get_full_state();
    auto ret = build_r_output<leapfrog::ChildModel, double>(full_state_output, save_steps);
    return ret;
}
