#pragma once

#include "../options.hpp"
#include "../generated/config_mixer.hpp"

namespace leapfrog {
namespace internal {

template<typename Config>
concept SpectrumPostHocCalculationsEnabled = RunDemographicProjection<Config> && RunHivSimulation<Config> && RunChildModel<Config> && RunSpectrumModel<Config>;

template<typename Config>
struct SpectrumPostHocCalculations {
  SpectrumPostHocCalculations(...) {};
};

template<SpectrumPostHocCalculationsEnabled Config>
struct SpectrumPostHocCalculations<Config> {
  using real_type = typename Config::real_type;
  using ModelVariant = typename Config::ModelVariant;
  using SS = Config::SS;
  using Pars = Config::Pars;
  using State = Config::State;
  using Intermediate = Config::Intermediate;
  using Args = Config::Args;

  // private members of this struct
  private:
  // state space
  static constexpr int NS = SS::NS;
  static constexpr int hDS = SS::hDS;
  static constexpr int hTS = SS::hTS;
  static constexpr int hAG = SS::hAG;
  static constexpr int pAG = SS::pAG;
  static constexpr int p_idx_hiv_first_adult = SS::p_idx_hiv_first_adult;

  // function args
  int t;
  const Pars& pars;
  const State& state_curr;
  State& state_next;
  Intermediate& intermediate;
  const Options<real_type>& opts;

  // only exposing the constructor and some methods
  public:
  SpectrumPostHocCalculations(Args& args):
    t(args.t),
    pars(args.pars),
    state_curr(args.state_curr),
    state_next(args.state_next),
    intermediate(args.intermediate),
    opts(args.opts)
  {};

  void run_spectrum_post_hoc_calulations() {
    calculate_nonaids_deaths();
  };

  // private methods that we don't want people to call
  private:
  void calculate_nonaids_deaths() {
    auto& n_ha = state_next.ha;
    auto& n_sp = state_next.sp;
    auto& i_sp = intermediate.sp;

    for (int g = 0; g < NS; ++g) {

      for (int ha = 0; ha < p_idx_hiv_first_adult; ++ha) {
        // Spectrum stores nonaids deaths by age but for children is always 0.
        // Write this out as 0 for children so we match Spectrum.
        n_sp.p_deaths_nonaids_artpop(ha, g) = 0.0;
        n_sp.p_deaths_nonaids_hivpop(ha, g) = 0.0;
      }

      for (int ha = p_idx_hiv_first_adult; ha < pAG; ++ha) {
        i_sp.hiv_art_adult_sa = 0.0;
        i_sp.hiv_untreated_adult_sa = 0.0;

        for (int hm = 0; hm < hDS; ++hm) {
          i_sp.hiv_untreated_adult_sa = i_sp.hiv_untreated_adult_sa + n_ha.h_hiv_adult(hm, ha - p_idx_hiv_first_adult, g);
          for (int hu = 0; hu < hTS; ++hu) {
            i_sp.hiv_art_adult_sa = i_sp.hiv_art_adult_sa + n_ha.h_art_adult(hu, hm, ha - p_idx_hiv_first_adult, g);
          }
        }

        if (i_sp.hiv_art_adult_sa + i_sp.hiv_untreated_adult_sa > 0) {
          i_sp.artcov_adult_sa = i_sp.hiv_art_adult_sa / (i_sp.hiv_art_adult_sa + i_sp.hiv_untreated_adult_sa);
        } else {
          i_sp.artcov_adult_sa = 0.0;
        }

        n_sp.p_deaths_nonaids_artpop(ha, g) = n_ha.p_hiv_pop_background_deaths(ha, g) * i_sp.artcov_adult_sa;
        n_sp.p_deaths_nonaids_hivpop(ha, g) = n_ha.p_hiv_pop_background_deaths(ha, g) * (1.0 - i_sp.artcov_adult_sa);
      }
    }
  };
};

}
}
