#pragma once

#include "state_space.hpp"

namespace leapfrog {

template<bool enable, typename Mixin>
struct Pair;

template<typename ...Ts>
struct SSMixer;

template<MV ModelVariant, typename Mixin, typename ...Ts>
struct SSMixer<ModelVariant, Pair<false, Mixin>, Ts...>: public SSMixer<ModelVariant, Ts...> {};

template<MV ModelVariant, typename Mixin, typename ...Ts>
struct SSMixer<ModelVariant, Pair<true, Mixin>, Ts...>: Mixin, public SSMixer<ModelVariant, Ts...> {};

template<MV ModelVariant>
struct SSMixer<ModelVariant>: public BaseSS {};

template<MV ModelVariant>
using SSMixed = SSMixer<
  ModelVariant,
  Pair<ModelVariant::run_demographic_projection, DpSS<ModelVariant>>,
  Pair<ModelVariant::run_hiv_simulation, HaSS<ModelVariant>>,
  Pair<ModelVariant::run_child_model, HcSS<ModelVariant>>
>;

}
