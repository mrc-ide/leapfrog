# leapfrog 0.0.3

* Implement Spectrum ART allocation.
  
  Spectrum allocates ART in a two step process: first, ART is allocated by CD4 category based
  on the 'expected mortality' and 'proportional to eligibility' weight. Second, within 
  CD4 categories, ART is allocated by age solely proportional to number in each age 
  group (propotional to eligibility).

* Patch ART dropout implementation. Spectrum converts input ART dropout percent to an 
  annual rate using [dropout rate] = -log(1.0 - [input percent]).

# leapfrog 0.0.2

* Add functions to read PMTCT and paediatric outputs from Spectrum .DP files.

# leapfrog 0.0.1

* Develop leapfrog model in templated C++ with R wrapper.

# leapfrog 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.
