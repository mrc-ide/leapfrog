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
  static constexpr auto hAG_span = SS::hAG_span;
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
    calculate_nonaids_excess_deaths();
  };

  // private methods that we don't want people to call
  private:
  void calculate_nonaids_deaths() {
    auto& n_ha = state_next.ha;
    auto& n_sp = state_next.sp;
    auto& i_sp = intermediate.sp;

    for (int g = 0; g < NS; ++g) {

      // Spectrum stores nonaids deaths by age but for children is always 0.
      // Only write values for a >= p_idx_hiv_first_adult
      // so that it is always 0 for children to match Spectrum.

      int a = p_idx_hiv_first_adult;
      for (int ha = 0; ha < hAG; ++ha) {
        i_sp.hiv_art_adult_sa = 0.0;
        i_sp.hiv_untreated_adult_sa = 0.0;

        for (int hm = 0; hm < hDS; ++hm) {
          i_sp.hiv_untreated_adult_sa += n_ha.h_hivpop(hm, ha, g);
          for (int hu = 0; hu < hTS; ++hu) {
            i_sp.hiv_art_adult_sa += n_ha.h_artpop(hu, hm, ha, g);
          }
        }

        if (i_sp.hiv_art_adult_sa + i_sp.hiv_untreated_adult_sa > 0) {
          i_sp.artcov_adult_sa = i_sp.hiv_art_adult_sa / (i_sp.hiv_art_adult_sa + i_sp.hiv_untreated_adult_sa);
        } else {
          i_sp.artcov_adult_sa = 0.0;
        }

        for (int i = 0; i < hAG_span[ha]; ++i, ++a) {
          n_sp.p_deaths_nonaids_artpop(a, g) = n_ha.p_deaths_background_hivpop(a, g) * i_sp.artcov_adult_sa;
          n_sp.p_deaths_nonaids_hivpop(a, g) = n_ha.p_deaths_background_hivpop(a, g) * (1.0 - i_sp.artcov_adult_sa);
        }
      }
    }
  };

  void calculate_nonaids_excess_deaths() {
    auto& n_ha = state_next.ha;
    auto& n_sp = state_next.sp;

    for (int g = 0; g < NS; ++g) {

      // Spectrum stores nonaids-excess deaths by age but for children
      // is always 0. Only write values for a >= p_idx_hiv_first_adult
      // so that it is always 0 for children to match Spectrum.

      int a = p_idx_hiv_first_adult;
      for (int ha = 0; ha < hAG; ++ha) {

        // Aggregate excess deaths and population in coarse age group
        auto excess_deaths_nonaids_no_art_ha = 0.0;
        auto excess_deaths_nonaids_on_art_ha = 0.0;
        auto hivpop_ha = 0.0;

        for (int hm = 0; hm < hDS; ++hm) {
          excess_deaths_nonaids_no_art_ha += n_ha.h_deaths_excess_nonaids_no_art(hm, ha, g);
          hivpop_ha += n_ha.h_hivpop(hm, ha, g);

          if (t > opts.ts_art_start) {
            for (int hu = 0; hu < hTS; ++hu) {
              excess_deaths_nonaids_on_art_ha += n_ha.h_deaths_excess_nonaids_on_art(hu, hm, ha, g);
              hivpop_ha += n_ha.h_artpop(hu, hm, ha, g);
            }
          }
        }


	      // Distribute deaths from coarse age group ha to single age group a proportional
	      // to distribution of HIV population in age group a
        for (int i = 0; i < hAG_span[ha]; ++i, ++a) {

          const auto hivpop_proportion_a = n_ha.p_hivpop(a, g) / hivpop_ha;
          n_sp.p_excess_deaths_nonaids_no_art(a, g) = excess_deaths_nonaids_no_art_ha * hivpop_proportion_a;

          if (t > opts.ts_art_start) {
            n_sp.p_excess_deaths_nonaids_on_art(a, g) += excess_deaths_nonaids_on_art_ha *  hivpop_proportion_a;
          }
        }
      }
    }
  }
};

}
}
