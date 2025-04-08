#pragma once

#include "cpp_adapter.hpp"

namespace leapfrog {
namespace internal {

template<typename ...Ts>
struct AdapterMixer;

template<typename real_type, MV ModelVariant>
struct AdapterMixer<real_type, ModelVariant> {
  using Config = ConfigMixer<real_type, ModelVariant>;

  static typename Config::Pars get_pars(const std::filesystem::path &input_dir, const Opts<real_type> &options, const int proj_years) {
    typename Config::Pars p = {}; return p;
  };

  static void write_output(std::filesystem::path& output_dir, const auto& state) {};
};

template<typename real_type, MV ModelVariant, typename Config, typename ...Ts>
struct AdapterMixer<real_type, ModelVariant, Pair<false, Config>, Ts...> : public AdapterMixer<real_type, ModelVariant, Ts...> {};

template<typename real_type1, MV ModelVariant1, typename ...Ts>
struct AdapterMixer<real_type1, ModelVariant1, Pair<true, DpConfig<real_type1, ModelVariant1>>, Ts...> {
  using real_type = real_type1;
  using ModelVariant = ModelVariant1;
  using CurrAdapter = DpAdapterCpp<real_type, ModelVariant>;
  using NextAdapterMixer = AdapterMixer<real_type, ModelVariant, Ts...>;
  using SS = SSMixed<ModelVariant>;
  using Config = ConfigMixer<real_type, ModelVariant, Pair<true, DpConfig<real_type1, ModelVariant1>>, Ts...>;

  static typename Config::Pars get_pars(const std::filesystem::path &input_dir, const Opts<real_type> &options, const int proj_years) {
    typename Config::Pars p = {
      NextAdapterMixer::get_pars(input_dir, options, proj_years),
      CurrAdapter::get_pars(input_dir, options, proj_years)
    };
    return p;
  };

  static void write_output(std::filesystem::path& output_dir, const typename Config::OutputState& state) {
    CurrAdapter::write_output(output_dir, state.dp);
    NextAdapterMixer::write_output(output_dir, state);
  };
};

template<typename real_type1, MV ModelVariant1, typename ...Ts>
struct AdapterMixer<real_type1, ModelVariant1, Pair<true, HaConfig<real_type1, ModelVariant1>>, Ts...> {
  using real_type = real_type1;
  using ModelVariant = ModelVariant1;
  using CurrAdapter = HaAdapterCpp<real_type, ModelVariant>;
  using NextAdapterMixer = AdapterMixer<real_type, ModelVariant, Ts...>;
  using SS = SSMixed<ModelVariant>;
  using Config = ConfigMixer<real_type, ModelVariant, Pair<true, HaConfig<real_type1, ModelVariant1>>, Ts...>;

  static typename Config::Pars get_pars(const std::filesystem::path &input_dir, const Opts<real_type> &options, const int proj_years) {
    typename Config::Pars p = {
      NextAdapterMixer::get_pars(input_dir, options, proj_years),
      CurrAdapter::get_pars(input_dir, options, proj_years)
    };
    return p;
  };

  static void write_output(std::filesystem::path& output_dir, const typename Config::OutputState& state) {
    CurrAdapter::write_output(output_dir, state.ha);
    NextAdapterMixer::write_output(output_dir, state);
  };
};

template<typename real_type1, MV ModelVariant1, typename ...Ts>
struct AdapterMixer<real_type1, ModelVariant1, Pair<true, HcConfig<real_type1, ModelVariant1>>, Ts...> {
  using real_type = real_type1;
  using ModelVariant = ModelVariant1;
  using CurrAdapter = HcAdapterCpp<real_type, ModelVariant>;
  using NextAdapterMixer = AdapterMixer<real_type, ModelVariant, Ts...>;
  using SS = SSMixed<ModelVariant>;
  using Config = ConfigMixer<real_type, ModelVariant, Pair<true, HcConfig<real_type1, ModelVariant1>>, Ts...>;

  static typename Config::Pars get_pars(const std::filesystem::path &input_dir, const Opts<real_type> &options, const int proj_years) {
    typename Config::Pars p = {
      NextAdapterMixer::get_pars(input_dir, options, proj_years),
      CurrAdapter::get_pars(input_dir, options, proj_years)
    };
    return p;
  };

  static void write_output(std::filesystem::path& output_dir, const typename Config::OutputState& state) {
    CurrAdapter::write_output(output_dir, state.hc);
    NextAdapterMixer::write_output(output_dir, state);
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
