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
* Refactor `OutputState` to take a struct of state-space dimensions instead of unpacking the subset of parameters we
  need. See https://github.com/mrc-ide/frogger/pull/12#discussion_r1245170775

## Leapfrog to Frogger glossary

Frogger naming rules

* `p_` prefix indicates this is age stratified by total single year i.e. same age stratification as total population
* `h_` prefix indicates this is age stratified by same age stratification as HIV adult population, this might be single
  year from
  15+ or could be coarse age groups
* `hc_` prefix indicates this is age stratified by the HIV child population
* `idx_` prefix means this is an index
* Uppercase means this is a dimension of the state space e.g. `NS` or `pAG`
* `ts_` prefix means this is a time step
* `hts_` prefix means this is an HIV time step (i.e. a time step of the inner HIV loop)

These prefixes can be merged e.g.

* `p_idx_hiv_first_adult` - Is the index of the first age group within singe-year ages to be considered an adult

### State space

In `StateSpace` struct

| Leapfrog | Frogger  | Details                                                        |
|----------|----------|----------------------------------------------------------------|
| NG       | NS       | Number of sexes                                                |
| pAG      | pAG      | Number of age groups in population, 81 for 0 to 80+            |
| hAG      | hAG      | Number of age groups in coarse stratified adult HIV population |
| hDS      | hDS      | Number of disease stages in adult HIV population               |
| hTS      | hTS      | Number of treatment stages in adult HIV population             |
| hAG_SPAN | hAG_span | Array of HIV age group sizes                                   |

#### Loop variable convention

* `NS` is `s`
* `pAG` is `a`
* `hAG` is `ha`
* `hDS` is `hd`
* `hTS` is `ht`

### Key indices in state space or other options controlling fit

In `Options` struct. Some of these come from data, some are static, some are an input option

| Leapfrog           | Frogger                | Details                                                |
|--------------------|------------------------|--------------------------------------------------------|
| hiv_steps_per_year | hts_per_year           | Number of HIV model time steps per year                |
| pIDX_FERT          | p_idx_fertility_first  | First age group index eligible for fertility           |
| pAG_FERT           | p_fertility_age_groups | Number of ages eligible for fertility                  |
| pIDX_HIVADULT      | p_idx_hiv_first_adult  | Index of the first age group to be considered an adult |
| t_ART_start        | ts_art_start           | Time step to start modelling ART                       |

### Input data

In 4 structs as part of the `Parameters`, named `Demography`, `Incidence`, `NaturalHistory` and `Art`

| Leapfrog           | Frogger                                  | Details                                                                                                                                                               |
|--------------------|------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| basepop            | demography.base_pop                      | Population data by age group and sex                                                                                                                                  |
| sx                 | demography.survival_probability          | Probability of surviving from age x to x+1 i.e. from 0 to 1, 1 to 2, ..., 79 to 80+ and 80+ to 80+                                                                    |
| netmigr            | demography.net_migration                 | Net migration by age group, sex and year                                                                                                                              |
| asfr               | demography.age_specific_fertility_rate   | Rate of live births in a year and the total female population of childbearing age                                                                                     |
| births_sex_prop    | demography.births_sex_prop               | Proportion of male and female births each year                                                                                                                        |
| incidinput         | incidence.rate                           | Incidence rate per year                                                                                                                                               |
| incrr_sex          | incidence.sex_rate_ratio                 | HIV incidence rate ratio for female : male age 15-49 years                                                                                                            |
| incrr_age          | incidence.age_rate_ratio                 | HIV incidence rate ratio by age (for each sex)                                                                                                                        |
| cd4_initdist       | natural_history.cd4_initial_distribution | Distribution of infections by CD4 category upon infection                                                                                                             |
| cd4_prog           | natural_history.cd4_progression          | Probability of progressing from 1 CD4 stage to the next by age and sex                                                                                                |
| cd4_mort           | natural_history.cd4_mortality            | Probability of mortality by CD4 stage, age and sex                                                                                                                    |
| scale_cd4_mort     | natural_history.scale_cd4_mortality      | If 1 then scale HIV related mortality (i.e. cd4_mortality) as a proportion of number of people with HIV and over the number with HIV and on ART at this disease stage |
| art_mort           | art.mortality_rate                       | Probability of mortality by treatment stage, CD4 stage, age and sex                                                                                                   |
| artmx_timerr       | art.mortaility_time_rate_ratio           | ART mortality rate ratio by year for <12 months and >12 months on ART                                                                                                 |
| art15plus_num      | art.adults_on_art                        | Time series of # or % of adult PLHIV on ART by sex                                                                                                                    |
| art15plus_isperc   | art.adults_on_art_is_percent             | Time series, TRUE if art15plus_num is a %, FALSE if it is a #                                                                                                         |
| artcd4elig_idx     | art.idx_hm_elig                          | The index of the CD4 count at which people are eligible for ART by time step                                                                                          |
| art_alloc_mxweight | art.initiation_mortality_weight          | Weighting for extent that expected mortality guides ART uptake                                                                                                        |
| art_dropout        | art.dropout                              | Annual ART dropout rate                                                                                                                                               |

### Outputs

Discussion: naming conventions

* Distinguish events vs. counts
* Distinguish stratification of array (single-year population pXX vs. HIV population hXX

| Leapfrog         | Frogger                    | Details                                                                  |
|------------------|----------------------------|--------------------------------------------------------------------------|
| totpop1          | p_total_pop                | Projected total population                                               |
| hivpop1          | p_hiv_pop                  | Projected HIV population                                                 |
| infections       | p_infections               | Projected number of new HIV infections by sex and age                    |
| hivstrat_adult   | h_hiv_adult                | Projected PLHIV not on ART by age, sex and CD4                           |
| artstrat_adult   | h_art_adult                | Projected PLHIV on ART by age, sex, CD4 and treatment stage              |
| births           | births                     | Projected number of births                                               |
| natdeaths        | p_total_pop_natural_deaths | Projected number of natural deaths                                       |
| natdeaths_hivpop | p_hiv_pop_natural_deaths   | Projected number of natural deaths within HIV population                 |
| hivdeaths        | p_hiv_deaths               | Projected HIV-related deaths by sex and age                              |
| aidsdeaths_noart | h_hiv_deaths_no_art        | Projected HIV-related deaths off ART by sex, age and CD4                 |
| aidsdeaths_art   | h_hiv_deaths_art           | Projected HIV-related deaths on ART by sex, age, CD4 and treatment stage |
| artinit          | h_art_initiation           | Projected ART initiations by sex, age and CD4                            |

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
| artelig_hahm              |                                 | Number of PLHIV eligible for ART by age and disease stage                 |
| Xart_15plus               |                                 |                                                                           |
| Xartelig_15plus           |                                 |                                                                           |
| expect_mort_artelig15plus |                                 | Total mortality amongst ART eligible population                           |
| artnum_hts                |                                 |                                                                           |
| artcov_hts                |                                 |                                                                           |
| curr_coverage             |                                 |                                                                           |
| artinit_hts               | art_initiation                  | Total number by sex to initiate on ART by sex                             |
| artinit_hahm              |                                 | Number who initiate treatment by age and disease stage                    |
| hivpop_ha(ha)             | hivpop_ha                       | HIV population by coarse HIV age group stratification                     |
| hivqx_ha                  |                                 |                                                                           |
| hivdeaths_a               |                                 |                                                                           |

## License

MIT © Imperial College of Science, Technology and Medicine
