from leapfrog_py import LeapfrogData, LeapfrogDataSingleYear, LeapfrogParameters, LeapfrogRange


def run_base_model(
    parameters: LeapfrogParameters,
    configuration: str,
    output_years: LeapfrogRange
) -> LeapfrogData: ...


def run_base_model_from_state(
    parameters: LeapfrogParameters,
    configuration: str,
    initial_state: LeapfrogDataSingleYear,
    simulation_start_year: int,
    output_years: LeapfrogRange
) -> LeapfrogData: ...


def run_base_model_single_year(
    parameters: LeapfrogParameters,
    configuration: str,
    initial_state: LeapfrogDataSingleYear,
    simulation_start_year: int
) -> LeapfrogDataSingleYear: ...
