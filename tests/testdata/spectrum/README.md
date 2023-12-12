# Spectrum test files

These are a series of test files constructed to 




Files were constructed using default data for Botswana
in Spectrum v5.13 on 12 February 2022.


## `bwa_aim-adult-no

This file includes 

ART is removed from the projecting

Paediatric HIV is removed by setting mother-to-child transmission probabilities to zero.

* Change all Adult ART percentages to zero.
* Change Adult ART "Percent lost to follow-up each year" to zero.
* Change all Child ART percentages to zero.
* Change child Percent receiving cotrimoxazole to zero.
* Set [Advanced Options] -> [MTCT transition probabilities] all equal zero.


## `v6.13/bwa_aim-adult-art-no-special-elig_v6.13_2022-04-18.PJNZ`

Paediatric HIV is removed by setting mother-to-child transmission probabilities to zero.

All ART inputs are percentages in the demo data file. Adult ART manually converted from 
percentages to counts for years 2000 to 2020 for consistency with actual usage.

The default file had no 'special population' ART eligibility. This is retained.

* Set [Advanced Options] -> [MTCT transition probabilities] all equal zero.
* Manually calculate number on ART 15+ and input for years 2000 to 2020.


## `v6.13/bwa2021_v6.13.pjnz`

* Created with Spectrum version 6.13

* Started with Botswana demo data for projection 1970 to 2030

* Adjusted adult ART eligibility:
  - 2008: change 200 to 250
  - 2012: change to <350 cells/mL
  - 2016: universal eligibility (999 cells/mL)
  - 2015: pregnant women eligible
  - 2013: co-infected eligible, 1.68%
  
* Adjusted paediatric ART eligibility
  - All children eligible: 2002: <12 months; 2012: <24 months; 2013: <60 months; 2017: <180 months
  
* Convert PMTCT inputs to numbers for up to year 2020; retain as percentages for year 2021 onwards
  - Use 'Convert values' button
  
* Convert child ART input to number age 0-14 through year 2019, age 0-4, 5-9, 10-14 for 2020 
  with distribution 15%, 30%, 55%, and percentage for 2021 onwards

* Convert adult ART input to number through year 2020, percentage for 2021 onwards


## `v6.28/` 

These files are copies of the `v6.13` test files re-projected and saved with Spectrum v6.28,
used for the 2023 UNAIDS Estimates.

A key update in Spectrum v6.28 is calendar-year demographic projection with net migration
all added at end-year, following UN WPP 2022 changes.
