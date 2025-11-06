from deepdiff import DeepDiff
from leapfrog_py import concat_on_time_dim, get_time_slice, run_model, run_model_from_state, run_model_single_year, read_h5_file


def assert_equal(obj1, obj2):
    assert len(DeepDiff(obj1, obj2, ignore_nan_inequality = True)) == 0


def test_adult_model_full_strat():
    parameters = read_h5_file("../r-package/tests/testthat/testdata/adult_parms_full.h5")
    ret = run_model(parameters)
    returned_vars = list(ret.keys())
    expected_vars = [
        "p_totpop", "births", "p_deaths_background_totpop", "p_hivpop",
        "p_deaths_background_hivpop", "h_hivpop", "h_artpop",
        "h_hiv_deaths_no_art", "p_infections", "h_hiv_deaths_art",
        "h_art_initiation", "h_deaths_excess_nonaids_no_art",
        "h_deaths_excess_nonaids_on_art", "p_deaths_excess_nonaids",
        "p_hiv_deaths", "p_net_migration_hivpop"
    ]
    returned_vars.sort()
    expected_vars.sort()
    assert_equal(returned_vars, expected_vars)


def test_child_model():
    parameters = read_h5_file("../r-package/tests/testthat/testdata/child_parms_full.h5")
    ret = run_model(parameters, "ChildModel")
    returned_vars = list(ret.keys())
    expected_vars = [
        "p_totpop", "births", "p_deaths_background_totpop", "p_hivpop",
        "p_deaths_background_hivpop", "h_hivpop", "h_artpop",
        "h_hiv_deaths_no_art", "p_infections", "h_hiv_deaths_art",
        "h_art_initiation", "h_deaths_excess_nonaids_no_art",
        "h_deaths_excess_nonaids_on_art", "p_deaths_excess_nonaids",
        "p_hiv_deaths", "p_net_migration_hivpop",
        "hiv_births", "hiv_births_by_mat_age",
        "hc1_hivpop", "hc2_hivpop",
        "hc1_artpop", "hc2_artpop",
        "hc1_noart_aids_deaths", "hc2_noart_aids_deaths",
        "hc1_art_aids_deaths", "hc2_art_aids_deaths",
        "hc_art_init", "hc_art_need_init", "ctx_need", "infection_by_type",
        "mtct_by_source_tr", "mtct_by_source_women",
        "mtct_by_source_hc_infections", "pmtct_coverage_at_delivery"
    ]
    returned_vars.sort()
    expected_vars.sort()
    assert_equal(returned_vars, expected_vars)


def test_child_model_running_twice_gives_same_result():
    parameters = read_h5_file("../r-package/tests/testthat/testdata/child_parms_full.h5")
    ret1 = run_model(parameters, "ChildModel")
    ret2 = run_model(parameters, "ChildModel")
    assert_equal(ret1, ret2)


def test_child_model_agrees_on_all_years_two_parts_single_year_runs():
    parameters = read_h5_file("../r-package/tests/testthat/testdata/child_parms_full.h5")
    ret_all_years = run_model(parameters, "ChildModel")

    ret_first_half = run_model(parameters, "ChildModel", range(1970, 2001))
    ret_second_half = run_model_from_state(parameters, "ChildModel", get_time_slice(ret_first_half, 30), 2000, range(2001, 2031))
    assert_equal(ret_all_years, concat_on_time_dim(ret_first_half, ret_second_half))

    ret_single_year = get_time_slice(run_model(parameters, "ChildModel", [1970]), 0)
    for year in range(1971, 2031):
        ret_single_year = run_model_single_year(parameters, "ChildModel", ret_single_year, year - 1)
        assert_equal(ret_single_year, get_time_slice(ret_all_years, year - 1970))
