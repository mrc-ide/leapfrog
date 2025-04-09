#pragma once

#include "../../options.hpp"
#include "cpp_adapter.hpp"
#include "../adapter_mixer.hpp"

namespace leapfrog {

template<typename real_type, internal::MV ModelVariant>
using AdapterCpp = internal::AdapterMixer<
  real_type, ModelVariant,
  internal::Triple<ModelVariant::run_demographic_projection, internal::DpConfig<real_type, ModelVariant>, internal::DpAdapterCpp<real_type, ModelVariant>>,
  internal::Triple<ModelVariant::run_hiv_simulation, internal::HaConfig<real_type, ModelVariant>, internal::HaAdapterCpp<real_type, ModelVariant>>,
  internal::Triple<ModelVariant::run_child_model, internal::HcConfig<real_type, ModelVariant>, internal::HcAdapterCpp<real_type, ModelVariant>>
>;

}
