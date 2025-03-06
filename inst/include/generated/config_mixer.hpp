#pragma once

#include "config.hpp"

namespace leapfrog {

template<typename ...Ts>
struct ConfigMixer;

template<typename real_type, MV ModelVariant>
struct ConfigMixer<real_type, ModelVariant> {
  struct Pars {};
  static Pars get_pars(const Rcpp::List &data, const Opts<real_type> &options, const int proj_years) {
    Pars p = {}; return p;
  };

  struct Intermediate {
    void reset() {};
  };

  struct State {
    State() {};
    State(const auto initial_state) {};
    void reset() {};
  };

  struct OutputState {
    OutputState(int output_years) {};
    void save_state(const size_t i, const auto &state) {};
  };

  static int get_build_output_size(int prev_size) {
    return prev_size;
  };
  static int build_output(Rcpp::List& ret, Rcpp::CharacterVector& names, int index, const auto& state, const size_t& output_years) {
    return index;
  };
};

template<typename real_type, MV ModelVariant, typename Config, typename ...Ts>
struct ConfigMixer<real_type, ModelVariant, Pair<false, Config>, Ts...> : public ConfigMixer<real_type, ModelVariant, Ts...> {};

#define ADD_CONFIG(VARIANT_NS, MEMBER_NS) \
template<typename real_type1, MV ModelVariant1, typename ...Ts> \
struct ConfigMixer<real_type1, ModelVariant1, Pair<true, VARIANT_NS ## Config<real_type1, ModelVariant1>>, Ts...> { \
  using real_type = real_type1; \
  using ModelVariant = ModelVariant1; \
  using CurrConfig = VARIANT_NS ## Config<real_type, ModelVariant>; \
  using NextConfigMixer = ConfigMixer<real_type, ModelVariant, Ts...>; \
  using SS = SSMixed<ModelVariant>; \
\
  struct Pars: public NextConfigMixer::Pars { \
    typename CurrConfig::Pars MEMBER_NS; \
  }; \
\
  static Pars get_pars(const Rcpp::List &data, const Opts<real_type> &options, const int proj_years) { \
    Pars p = { \
      NextConfigMixer::get_pars(data, options, proj_years), \
      CurrConfig::get_pars(data, options, proj_years) \
    }; \
    return p; \
  }; \
\
  struct Intermediate: public NextConfigMixer::Intermediate { \
    typename CurrConfig::Intermediate MEMBER_NS; \
\
    Intermediate(): \
      NextConfigMixer::Intermediate(), \
      MEMBER_NS() {}; \
\
    void reset() { \
      NextConfigMixer::Intermediate::reset(); \
      MEMBER_NS.reset(); \
    }; \
  }; \
\
  struct State: public NextConfigMixer::State { \
    typename CurrConfig::State MEMBER_NS; \
\
    State(): \
      NextConfigMixer::State(), \
      MEMBER_NS() {}; \
\
    State(const typename CurrConfig::State& initial_state): \
      NextConfigMixer::State(initial_state), \
      MEMBER_NS(initial_state) {}; \
\
    void reset() { \
      NextConfigMixer::State::reset(); \
      MEMBER_NS.reset(); \
    }; \
  }; \
\
  struct OutputState: public NextConfigMixer::OutputState { \
    typename CurrConfig::OutputState MEMBER_NS; \
\
    OutputState(int output_years): \
      NextConfigMixer::OutputState(output_years), \
      MEMBER_NS(output_years) {}; \
\
    void save_state(const size_t i, const State &state) { \
      NextConfigMixer::OutputState::save_state(i, state); \
      MEMBER_NS.save_state(i, state.MEMBER_NS); \
    }; \
  }; \
\
  static int get_build_output_size(int prev_size) { \
    int curr_size = NextConfigMixer::get_build_output_size(prev_size); \
    return CurrConfig::get_build_output_size(curr_size); \
  }; \
\
  static int build_output(Rcpp::List& ret, Rcpp::CharacterVector& names, int index, const OutputState& state, const size_t& output_years) { \
    int new_index = CurrConfig::build_output(ret, names, index, state.MEMBER_NS, output_years); \
    return NextConfigMixer::build_output(ret, names, new_index, state, output_years); \
  }; \
\
  using Options = Opts<real_type>; \
  struct Args { \
    int t; \
    const Pars& pars; \
    const State& state_curr; \
    State& state_next; \
    Intermediate& intermediate; \
    const Options& opts; \
  }; \
};

ADD_CONFIG(Dp, dp)
ADD_CONFIG(Ha, ha)
ADD_CONFIG(Hc, hc)

template<typename real_type, MV ModelVariant>
using ConfigMixed = ConfigMixer<
    real_type, ModelVariant,
    Pair<ModelVariant::run_demographic_projection, DpConfig<real_type, ModelVariant>>,
    Pair<ModelVariant::run_hiv_simulation, HaConfig<real_type, ModelVariant>>,
    Pair<ModelVariant::run_child_model, HcConfig<real_type, ModelVariant>>
>;

}
