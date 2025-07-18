def run_base_model(
    parameters: dict,
    configuration: str,
    output_years: list
) -> dict: ...

def run_base_model_from_state(
    parameters: dict,
    configuration: str,
    initial_state: dict,
    simulation_start_year: int,
    output_years: list
) -> dict: ...

def run_base_model_single_year(
    parameters: dict,
    configuration: str,
    initial_state: dict,
    simulation_start_year: int
) -> dict: ...
