#include <Rcpp.h>

#include "leapfrog-raw.h"

//' Simulate leapfrog model
//'
//' @param demp list of demographic input parameters (TODO: document)
//' @param projp list of HIV projection parameters (TODO: document)
//' @param hiv_strat stratification of HIV population, either "full"
//'   (default; single-year ages) or "coarse" (aggregated age groups). 
//' @param hiv_steps_per_year number of Euler integration steps per year
//'   for HIV progression; default 10.
//'
//' @details
//' The first year of `sx`, `asfr`, `srb`, and `netmig` is not used. This is assumed
//' to apply to the base year population (consistent with Spectrum).
//'
//' @export
//' 
// [[Rcpp::export]]
Rcpp::List
leapfrogR(const Rcpp::List& demp,
	  const Rcpp::List& projp,
	  const Rcpp::String hiv_strat = "full",
	  const int hiv_steps_per_year = 10) {

  using namespace Rcpp;
  
  NumericVector Sx = demp["Sx"];
  Dimension d = Sx.attr("dim");
  const size_t proj_years = d[2];
  const int NG = 2;
  const int pAG = 81;
  const int pIDX_HIVADULT = 15;
  const int pIDX_FERT = 15;
  const int pAG_FERT = 35;
  const int hAG_COARSE = 9;
  const int hAG_FULL = 66;
  const int hDS = 7;
  const int hTS = 3;

  int hAG;
  if (hiv_strat == "full") {
    hAG = hAG_FULL;
  } else if (hiv_strat == "coarse") {
    hAG = hAG_COARSE;
  } else {
    Rf_error("hiv_strat \"%s\" not found. Please select \"full\" or \"coarse\".\n", hiv_strat.get_cstring());
  }
    
  // allocate memory for return object
  NumericVector totpop1(pAG * NG * proj_years);
  totpop1.attr("dim") = NumericVector::create(pAG, NG, proj_years);

  NumericVector hivpop1(pAG * NG * proj_years);
  hivpop1.attr("dim") = NumericVector::create(pAG, NG, proj_years);

  NumericVector infections(pAG * NG * proj_years);
  infections.attr("dim") = NumericVector::create(pAG, NG, proj_years);

  NumericVector hivstrat_adult(hDS * hAG * NG * proj_years);
  hivstrat_adult.attr("dim") = NumericVector::create(hDS, hAG, NG, proj_years);

  NumericVector artstrat_adult(hTS * hDS * hAG * NG * proj_years);
  artstrat_adult.attr("dim") = NumericVector::create(hTS, hDS, hAG, NG, proj_years);

  NumericVector natdeaths(pAG * NG * proj_years);
  natdeaths.attr("dim") = NumericVector::create(pAG, NG, proj_years);

  NumericVector natdeaths_hivpop(pAG * NG * proj_years);
  natdeaths_hivpop.attr("dim") = NumericVector::create(pAG, NG, proj_years);

  NumericVector hivdeaths(pAG * NG * proj_years);
  hivdeaths.attr("dim") = NumericVector::create(pAG, NG, proj_years);

  NumericVector aidsdeaths_noart(hDS * hAG * NG * proj_years);
  aidsdeaths_noart.attr("dim") = NumericVector::create(hDS, hAG, NG, proj_years);

  NumericVector aidsdeaths_art(hTS * hDS * hAG * NG * proj_years);
  aidsdeaths_art.attr("dim") = NumericVector::create(hTS, hDS, hAG, NG, proj_years);

  NumericVector artinit(hDS * hAG * NG * proj_years);
  artinit.attr("dim") = NumericVector::create(hDS, hAG, NG, proj_years);
 
  if (hAG == hAG_FULL) {
    leapfrog_sim<double, NG, pAG, pIDX_FERT, pAG_FERT,
		 pIDX_HIVADULT, hAG_FULL, hDS, hTS>
      (REAL(demp["basepop"]),
       REAL(demp["Sx"]),
       REAL(demp["netmigr_adj"]),
       REAL(demp["asfr"]),
       REAL(demp["births_sex_prop"]),
       REAL(projp["incidinput"]),
       REAL(projp["incrr_sex"]),
       REAL(projp["incrr_age"]),
       REAL(projp["cd4_initdist_full"]),
       REAL(projp["cd4_prog_full"]),
       REAL(projp["cd4_mort_full"]),
       REAL(projp["art_mort_full"]),
       REAL(projp["artmx_timerr"]),
       REAL(projp["art15plus_num"]),
       LOGICAL(projp["art15plus_isperc"]),
       INTEGER(projp["artcd4elig_idx"]),
       *INTEGER(projp["art_alloc_method"]),
       *REAL(projp["art_alloc_mxweight"]),
       *INTEGER(projp["scale_cd4_mort"]),
       REAL(projp["art_dropout"]),
       proj_years,
       hiv_steps_per_year,
       *INTEGER(projp["t_ART_start"]) - 1, // 0-based indexing vs. R 1-based
       INTEGER(projp["hAG_SPAN_full"]),
       REAL(totpop1),
       REAL(hivpop1),
       REAL(infections),
       REAL(hivstrat_adult),
       REAL(artstrat_adult),
       REAL(natdeaths),
       REAL(natdeaths_hivpop),
       REAL(hivdeaths),
       REAL(aidsdeaths_noart),
       REAL(aidsdeaths_art),
       REAL(artinit));
  } else if (hAG == hAG_COARSE) {
    leapfrog_sim<double, NG, pAG, pIDX_FERT, pAG_FERT,
		 pIDX_HIVADULT, hAG_COARSE, hDS, hTS>
      (REAL(demp["basepop"]),
       REAL(demp["Sx"]),
       REAL(demp["netmigr_adj"]),
       REAL(demp["asfr"]),
       REAL(demp["births_sex_prop"]),
       REAL(projp["incidinput"]),
       REAL(projp["incrr_sex"]),
       REAL(projp["incrr_age"]),
       REAL(projp["cd4_initdist_coarse"]),
       REAL(projp["cd4_prog_coarse"]),
       REAL(projp["cd4_mort_coarse"]),
       REAL(projp["art_mort_coarse"]),
       REAL(projp["artmx_timerr"]),
       REAL(projp["art15plus_num"]),
       LOGICAL(projp["art15plus_isperc"]),
       INTEGER(projp["artcd4elig_idx"]),
       *INTEGER(projp["art_alloc_method"]),
       *REAL(projp["art_alloc_mxweight"]),
       *INTEGER(projp["scale_cd4_mort"]),
       REAL(projp["art_dropout"]),       
       proj_years,
       hiv_steps_per_year,
       *INTEGER(projp["t_ART_start"]) - 1,  // 0-based indexing vs. R 1-based
       INTEGER(projp["hAG_SPAN_coarse"]),
       REAL(totpop1),
       REAL(hivpop1),
       REAL(infections),
       REAL(hivstrat_adult),
       REAL(artstrat_adult),
       REAL(natdeaths),
       REAL(natdeaths_hivpop),
       REAL(hivdeaths),
       REAL(aidsdeaths_noart),
       REAL(aidsdeaths_art),
       REAL(artinit));
  } else {
    Rf_error("Invalid HIV stratification age groups (hAG)");
  }

  List ret = List::create(_("totpop1") = totpop1,
			  _("hivpop1") = hivpop1,
			  _("hivstrat_adult") = hivstrat_adult,
			  _("artstrat_adult") = artstrat_adult,
			  _("infections") = infections,
			  _("natdeaths") = natdeaths,
			  _("natdeaths_hivpop") = natdeaths_hivpop,
			  _("hivdeaths") = hivdeaths,
			  _("aidsdeaths_noart") = aidsdeaths_noart,
			  _("aidsdeaths_art") = aidsdeaths_art,
			  _("artinit") = artinit);
				      
  return ret;
}
