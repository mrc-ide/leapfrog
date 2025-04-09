#pragma once

#include "../../options.hpp"
#include "r_adapter.hpp"
#include "../adapter_mixer.hpp"

namespace leapfrog {

template<typename real_type, internal::MV ModelVariant>
using AdapterR = internal::AdapterMixer<
  real_type, ModelVariant,
  internal::Triple<ModelVariant::run_demographic_projection, internal::DpConfig<real_type, ModelVariant>, internal::DpAdapterR<real_type, ModelVariant>>,
  internal::Triple<ModelVariant::run_hiv_simulation, internal::HaConfig<real_type, ModelVariant>, internal::HaAdapterR<real_type, ModelVariant>>,
  internal::Triple<ModelVariant::run_child_model, internal::HcConfig<real_type, ModelVariant>, internal::HcAdapterR<real_type, ModelVariant>>
>;

}
