#pragma once

#include "r_adapter.hpp"

namespace leapfrog {

namespace internal {

template<typename ...Ts>
struct AdapterMixer;

template<typename real_type, MV ModelVariant>
struct AdapterMixer<real_type, ModelVariant> {
  using Config = ConfigMixer<real_type, ModelVariant>;

  static typename Config::Pars get_pars(const Rcpp::List &data, const Opts<real_type> &options, const int proj_years) {
    typename Config::Pars p = {}; return p;
  };

  static int build_output(Rcpp::List& ret, Rcpp::CharacterVector& names, int index, const auto& state, const size_t& output_years) {
    return index;
  };
};

template<typename real_type, MV ModelVariant, typename Config, typename ...Ts>
struct AdapterMixer<real_type, ModelVariant, Pair<false, Config>, Ts...> : public AdapterMixer<real_type, ModelVariant, Ts...> {};

template<typename real_type1, MV ModelVariant1, typename ...Ts>
struct AdapterMixer<real_type1, ModelVariant1, Pair<true, DpConfig<real_type1, ModelVariant1>>, Ts...> {
  using real_type = real_type1;
  using ModelVariant = ModelVariant1;
  using CurrAdapter = DpAdapterR<real_type, ModelVariant>;
  using NextAdapterMixer = AdapterMixer<real_type, ModelVariant, Ts...>;
  using SS = SSMixed<ModelVariant>;
  using Config = ConfigMixer<real_type, ModelVariant, Pair<true, DpConfig<real_type1, ModelVariant1>>, Ts...>;

  static typename Config::Pars get_pars(const Rcpp::List &data, const Opts<real_type> &options, const int proj_years) {
    typename Config::Pars p = {
      NextAdapterMixer::get_pars(data, options, proj_years),
      CurrAdapter::get_pars(data, options, proj_years)
    };
    return p;
  };

  static int build_output(Rcpp::List& ret, Rcpp::CharacterVector& names, int index, const typename Config::OutputState& state, const size_t& output_years) {
    int new_index = CurrAdapter::build_output(ret, names, index, state.dp, output_years);
    return NextAdapterMixer::build_output(ret, names, new_index, state, output_years);
  };
};

template<typename real_type1, MV ModelVariant1, typename ...Ts>
struct AdapterMixer<real_type1, ModelVariant1, Pair<true, HaConfig<real_type1, ModelVariant1>>, Ts...> {
  using real_type = real_type1;
  using ModelVariant = ModelVariant1;
  using CurrAdapter = HaAdapterR<real_type, ModelVariant>;
  using NextAdapterMixer = AdapterMixer<real_type, ModelVariant, Ts...>;
  using SS = SSMixed<ModelVariant>;
  using Config = ConfigMixer<real_type, ModelVariant, Pair<true, HaConfig<real_type1, ModelVariant1>>, Ts...>;

  static typename Config::Pars get_pars(const Rcpp::List &data, const Opts<real_type> &options, const int proj_years) {
    typename Config::Pars p = {
      NextAdapterMixer::get_pars(data, options, proj_years),
      CurrAdapter::get_pars(data, options, proj_years)
    };
    return p;
  };

  static int build_output(Rcpp::List& ret, Rcpp::CharacterVector& names, int index, const typename Config::OutputState& state, const size_t& output_years) {
    int new_index = CurrAdapter::build_output(ret, names, index, state.ha, output_years);
    return NextAdapterMixer::build_output(ret, names, new_index, state, output_years);
  };
};

template<typename real_type1, MV ModelVariant1, typename ...Ts>
struct AdapterMixer<real_type1, ModelVariant1, Pair<true, HcConfig<real_type1, ModelVariant1>>, Ts...> {
  using real_type = real_type1;
  using ModelVariant = ModelVariant1;
  using CurrAdapter = HcAdapterR<real_type, ModelVariant>;
  using NextAdapterMixer = AdapterMixer<real_type, ModelVariant, Ts...>;
  using SS = SSMixed<ModelVariant>;
  using Config = ConfigMixer<real_type, ModelVariant, Pair<true, HcConfig<real_type1, ModelVariant1>>, Ts...>;

  static typename Config::Pars get_pars(const Rcpp::List &data, const Opts<real_type> &options, const int proj_years) {
    typename Config::Pars p = {
      NextAdapterMixer::get_pars(data, options, proj_years),
      CurrAdapter::get_pars(data, options, proj_years)
    };
    return p;
  };

  static int build_output(Rcpp::List& ret, Rcpp::CharacterVector& names, int index, const typename Config::OutputState& state, const size_t& output_years) {
    int new_index = CurrAdapter::build_output(ret, names, index, state.hc, output_years);
    return NextAdapterMixer::build_output(ret, names, new_index, state, output_years);
  };

};

}

template<typename real_type, internal::MV ModelVariant>
using Adapter = internal::AdapterMixer<
  real_type, ModelVariant,
  internal::Pair<ModelVariant::run_demographic_projection, internal::DpConfig<real_type, ModelVariant>>,
  internal::Pair<ModelVariant::run_hiv_simulation, internal::HaConfig<real_type, ModelVariant>>,
  internal::Pair<ModelVariant::run_child_model, internal::HcConfig<real_type, ModelVariant>>
>;

}
