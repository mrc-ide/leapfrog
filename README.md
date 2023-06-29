# frogger

<!-- badges: start -->
[![Project Status: Concept – Minimal or no implementation has been done yet, or the repository is only intended to be a limited example, demo, or proof-of-concept.](https://www.repostatus.org/badges/latest/concept.svg)](https://www.repostatus.org/#concept)
[![R build status](https://github.com/mrc-ide/frogger/workflows/R-CMD-check/badge.svg)](https://github.com/mrc-ide/frogger/actions)
[![codecov.io](https://codecov.io/github/mrc-ide/frogger/coverage.svg?branch=main)](https://codecov.io/github/mrc-ide/frogger?branch=main)
[![CodeFactor](https://www.codefactor.io/repository/github/mrc-ide/frogger/badge)](https://www.codefactor.io/repository/github/mrc-ide/frogger)
<!-- badges: end -->

## Installation

To install `frogger`:

```r
remotes::install_github("mrc-ide/frogger", upgrade = FALSE)
```

## TODO

* Improve variable names!
* Restructuring the model code to identify more common code
    * There are examples like general demographic projection and hiv population demographic projection which are running
      similar processes like ageing, non HIV mortality, migration. We should be able to write a function for e.g. ageing
      which we can run on each of our population matrices. Even for the HIV and ART stratified we can add overloaded
      function to work with higher dimension data
* Add a test that checks that no `double`s are used in `inst/include` dir. We should be using templated `real_type` for
  TMB
* Add a broad level overview of the algorithm - is there a diagram available?
* Add a test which proves that C++ code can be compiled standalone, one exists in dust that can be used as inspriration
* Tidy data copying a tensor2 to matrix and a tensor n to array utility functions which should be used in C++ R
  interface "src/frogger.cpp"
* Ensure input data is not copied, we can read from the data that R owns (as long as don't write to it)
* Convert string flags to the model to enums if we need to switch on them in several places. This should make it easier
  to reason about in C++ world and isolate the string checking to a single place
* Update `gender` terminology to sex
* Previously `hiv_negative_pop` was fixed size by having dimensions specified by template, how much does this speed up
  the code? Is there a better way to do this?
* Tidy up confusing looping see https://github.com/mrc-ide/frogger/pull/7#discussion_r1217847753
* Add R casting helpers which return better errors than Rcpp
  see https://github.com/mrc-ide/frogger/pull/7#discussion_r1217884684
* Add a helper to do 0 to base 1 conversion and check upper bounds
  see https://github.com/mrc-ide/frogger/pull/7#discussion_r1217888684
* Review what we pass as parameters - can some of these be computed in the struct ctor?
  e.g. https://github.com/mrc-ide/frogger/pull/7#discussion_r1217890466
* Add function for calculating incidence rate per sex https://github.com/mrc-ide/frogger/pull/7#discussion_r1217848792
* Make fit work with coarse ages (at the moment not reading all of the coarse stratified data)
* Remove duplicate reading of `hAG_SPAN_full`, read this as `hiv_age_groups_span` and `age_groups_hiv_span`

## Leapfrog to Frogger glossary

### State space

In `StateSpace` struct

| Leapfrog           | Frogger            | Details                                             |
|--------------------|--------------------|-----------------------------------------------------|
| NG                 | NG                 | Number of genders                                   |
| pAG                | pAG                | Number of age groups in population, 81 for 0 to 80+ |
| hAG                | hAG                | Number of age groups in HIV population              |
| hDS                | hDS                | Number of disease stages                            |
| hTS                | hTS                | Number of treatment stages                          |
| sim_years          | sim_years          | Number of simulation years to run                   |
| hiv_steps_per_year | hiv_steps_per_year | Number of HIV model time steps per year             |

### Key indices in state space

In `?` struct

| Leapfrog      | Frogger                   | Details                                                |
|---------------|---------------------------|--------------------------------------------------------|
| pIDX_FERT     | fertility_first_age_group | First age group index eligible for fertility           |
| pAG_FERT      | age_groups_fertiity       | Number of ages eligible for fertility                  |
| pIDX_HIVADULT | hiv_adult_first_age_group | Index of the first age group to be considered an adult |
| t_ART_start   | time_art_start            | Time step to start modelling ART                       |
| hAG_SPAN      | hiv_age_groups_span       | Array of HIV age group sizes                           |

### Input data

In 4 structs as part of the `Parameters`, named `Demography`, `Incidence`, `NaturalHistory` and `Art`

| Leapfrog           | Frogger                                  | Details                                                                                                                                                               |
|--------------------|------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| basepop            | demography.base_pop                      | Population data by age group and sex                                                                                                                                  |
| sx                 | demography.survival                      | Probability of surviving between ages, from 0 to 1, 1 to 2, ..., 79 to 80+ and 80+ to 80+                                                                             |
| netmigr            | demography.net_migration                 | Net migration by age group, sex and year                                                                                                                              |
| asfr               | demography.fertility_ratio               | Ratio of number of live births in a year and the whole female population of childbearing age                                                                          |
| births_sex_prop    | demography.births_sex_prop               | Proportion of male and female births each year                                                                                                                        |
| incidinput         | incidence.rate                           | Incidence rate per year                                                                                                                                               |
| incrr_sex          | incidence.sex_rate_ratio                 | HIV incidence rate ratio for female : male age 15-49 years                                                                                                            |
| incrr_age          | incidence.age_rate_ratio                 | HIV incidence rate ratio by age (for each sex)                                                                                                                        |
| cd4_initdist       | natural_history.cd4_initial_distribution | Distribution of infections by CD4 category upon infection                                                                                                             |
| cd4_prog           | natural_history.cd4_progression          | Probability of progressing from 1 CD4 stage to the next by age and sex                                                                                                |
| cd4_mort           | natural_history.cd4_mortality            | Probability of mortality by CD4 stage, age and sex                                                                                                                    |
| scale_cd4_mort     | natural_history.scale_cd4_mortality      | If 1 then scale HIV related mortality (i.e. cd4_mortality) as a proportion of number of people with HIV and over the number with HIV and on ART at this disease stage |
| art_mort           | art.mortality                            | Probability of mortality by treatment stage, CD4 stage, age and sex                                                                                                   |
| artmx_timerr       | ?                                        |                                                                                                                                                                       |
| art15plus_num      | ?                                        |                                                                                                                                                                       |
| art15plus_isperc   | ?                                        | TRUE if art                                                                                                                                                           |
| artcd4elig_idx     | ?                                        | The index of the CD4 count at which people are eligible for ART by time step                                                                                          |
| art_alloc_mxweight | ?                                        |                                                                                                                                                                       |
| art_dropout        | art.dropout                              | ART dropout rate at each time step                                                                                                                                    |

### Outputs

_Discussion: naming conventions

* Distinguish events vs. counts
* Distinguish stratification of array (single-year population pXX vs. HIV population hXX

| Leapfrog         | Frogger              | Details                                                  |
|------------------|----------------------|----------------------------------------------------------|
| totpop1          | p_total_population   | Projected total population                               |
| hivpop1          | p_hiv_population     | Projected HIV population                                 |
| infections       | p_infections         |                                                          |
| hivstrat_adult   | h_hiv_strat_adult    |                                                          |
| artstrat_adult   | h_art_strat_adult    |                                                          |
| births           | births               | Projected number of births                               |
| natdeaths        | p_natural_deaths     | Projected number of natural deaths                       |
| natdeaths_hivpop | p_natural_deaths_hiv | Projected number of natural deaths within HIV population |
| hivdeaths        | p_hiv_deaths         |                                                          |
| aidsdeaths_noart | h_aids_deaths_no_art |                                                          |
| aidsdeaths_art   | h_aids_deaths_art    |                                                          |
| artinit          | h_art_initiation     |                                                          |

### Internal

| Leapfrog                  | Frogger                         | Details                                                                   |
|---------------------------|---------------------------------|---------------------------------------------------------------------------|
| ART0MOS                   | ?                               |                                                                           |
| pIDX_INCIDPOP             | adult_incidence_first_age_group | Index of youngest age that is reflected in the adult incidence input      |
| pAG_INCIDPOP              | ?                               |                                                                           |
| hAG_15PLUS                | age_groups_hiv_15plus           |                                                                           |
| hIDX_15PLUS               | hIDX_15PLUS                     |                                                                           |
| h_art_stage_dur           | ?                               |                                                                           |
| everARTelig_idx           | everARTelig_idx                 |                                                                           |
| migrate_ag                | migration_rate                  | Rate people migrate into population by age and sex                        |
| sx_netmig                 | survival_netmig                 |                                                                           |
| births_sex                | births_sex                      |                                                                           |
| migrate_a0                | migration_rate_a0               |                                                                           |
| hiv_ag_prob               | hiv_age_up_prob                 | Probability of aging from one group to the next in HIV age stratification |
| hivpop_ha(ha, g)          | hiv_population_coarse_ages      |                                                                           |
| netmig_ag                 | hiv_net_migration               |                                                                           |
| deathsmig_ha              | deaths_migrate                  |                                                                           |
| deathmigrate_ha           | deaths_migrate_rate             |                                                                           |
| cd4elig_idx               | ?                               | Index of the maximum CD4 count category eligible for ART                  |
| anyelig_idx               | anyelig_idx                     |                                                                           |
| everARTelig_idx           | ?                               |                                                                           |
| infections_ts             | infections_ts                   | Infections occurring at a specific time step                              |
| hivn_ag                   | hiv_negative_pop                | HIV negative population for a given age and sex                           |
| Xhivn                     | hiv_neg_aggregate               |                                                                           |
| Xhivn_incagerr            | ?                               |                                                                           |
| incrate_i                 |                                 | Incidence rate input                                                      |
| incrate_g                 | incidence_rate_sex              | Incidence rate by sex                                                     |
| hivdeaths_ha              | hiv_deaths_age_sex              |                                                                           |
| grad                      | grad                            | Movement between cd4 categories at a given time step                      |
| artpop_hahm               |                                 |                                                                           |
| cd4mx_scale               |                                 |                                                                           |
| infections_a              | infections_a                    |                                                                           |
| infections_ha             | infections_ha                   |                                                                           |
| gradART                   |                                 |                                                                           |
| artelig_hahm              |                                 |                                                                           |
| Xart_15plus               |                                 |                                                                           |
| Xartelig_15plus           |                                 |                                                                           |
| expect_mort_artelig15plus |                                 |                                                                           |
| artnum_hts                |                                 |                                                                           |
| artcov_hts                |                                 |                                                                           |
| curr_coverage             |                                 |                                                                           |
| artinit_hts               |                                 |                                                                           |
| artinit_hahm              |                                 |                                                                           |
| hivpop_ha(ha)             | hivpop_ha                       | HIV population by coarse HIV age group stratification                     |
| hivqx_ha                  |                                 |                                                                           |
| hivdeaths_a               |                                 |                                                                           |

## License

MIT © Imperial College of Science, Technology and Medicine
