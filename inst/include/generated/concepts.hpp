#pragma once

namespace leapfrog {

template<typename ModelVariant>
concept MV = requires {
    { ModelVariant::run_demographic_projection } -> std::convertible_to<bool>;
    { ModelVariant::run_hiv_simulation } -> std::convertible_to<bool>;
    { ModelVariant::use_full_stratification } -> std::convertible_to<bool>;
    { ModelVariant::run_child_model } -> std::convertible_to<bool>;
};

template<typename Config>
concept RunDemographicProjection = MV<typename Config::ModelVariant> && Config::ModelVariant::run_demographic_projection;

template<typename Config>
concept RunHivSimulation = MV<typename Config::ModelVariant> && Config::ModelVariant::run_hiv_simulation;

template<typename Config>
concept UseFullStratification = MV<typename Config::ModelVariant> && Config::ModelVariant::use_full_stratification;

template<typename Config>
concept RunChildModel = MV<typename Config::ModelVariant> && Config::ModelVariant::run_child_model;

template<typename Config>
concept GeneralDemographicProjectionEnabled = RunDemographicProjection<Config>;

template<typename Config>
concept HivDemographicProjectionEnabled = RunDemographicProjection<Config> && RunHivSimulation<Config>;

template<typename Config>
concept HivModelSimulationEnabled = RunDemographicProjection<Config> && RunHivSimulation<Config>;

template<typename Config>
concept ChildModelSimulationEnabled = RunDemographicProjection<Config> && RunHivSimulation<Config> && RunChildModel<Config>;

}
