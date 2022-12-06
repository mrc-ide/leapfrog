#include <Rcpp.h>

#include "leapfrog-raw.h"
#include "ListBuilder.h" // alternative for construct lists longer than 20 elements


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
  const int hDS_adol = 6;
  const int trans = 4;
  const int tx_time = 3;
  const int hTS = 3;
  const int ctx_effect = 0.33;

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

  NumericVector hivnpop1(pAG * NG * proj_years);
  hivnpop1.attr("dim") = NumericVector::create(pAG, NG, proj_years);

  NumericVector infections(pAG * NG * proj_years);
  infections.attr("dim") = NumericVector::create(pAG, NG, proj_years);

  NumericVector hivstrat_adult(hDS * hAG * NG * proj_years);
  hivstrat_adult.attr("dim") = NumericVector::create(hDS, hAG, NG, proj_years);

  NumericVector artstrat_adult(hTS * hDS * hAG * NG * proj_years);
  artstrat_adult.attr("dim") = NumericVector::create(hTS, hDS, hAG, NG, proj_years);
  
  NumericVector hivstrat_paeds(trans * hDS * pIDX_HIVADULT * NG * proj_years);
  hivstrat_paeds.attr("dim") = NumericVector::create(hDS, trans, pIDX_HIVADULT, NG, proj_years);
  
  NumericVector artstrat_paeds(hTS * hDS * pIDX_HIVADULT * NG * proj_years);
  artstrat_paeds.attr("dim") = NumericVector::create(hTS, hDS, pIDX_HIVADULT, NG, proj_years);
  
  NumericVector artelig_paeds(hDS * trans * pIDX_HIVADULT * NG * proj_years);
  artelig_paeds.attr("dim") = NumericVector::create(hDS, trans, pIDX_HIVADULT, NG, proj_years);
  
  NumericVector coarse_totpop1(hAG * NG * proj_years);
  coarse_totpop1.attr("dim") = NumericVector::create(hAG, NG, proj_years);

  NumericVector births(proj_years);
  
  NumericVector hiv_births(proj_years);
  

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
  
  NumericVector aidsdeaths_noart_paed(hDS * pIDX_HIVADULT * NG * proj_years);
  aidsdeaths_noart_paed.attr("dim") = NumericVector::create(hDS, pIDX_HIVADULT, NG, proj_years);
  
  NumericVector aidsdeaths_art_paed(hTS * hDS * pIDX_HIVADULT * NG * proj_years);
  aidsdeaths_art_paed.attr("dim") = NumericVector::create(hTS, hDS, pIDX_HIVADULT, NG, proj_years);

  NumericVector artinit(hDS * hAG * NG * proj_years);
  artinit.attr("dim") = NumericVector::create(hDS, hAG, NG, proj_years);
  
  
  NumericVector artnum_paed(proj_years);
  
  NumericVector deaths_paeds(trans * hDS * pIDX_HIVADULT * NG * proj_years);
  deaths_paeds.attr("dim") = NumericVector::create(hDS, trans,  pIDX_HIVADULT, NG, proj_years);
  
  NumericVector deaths_paeds_art(trans * hDS * pIDX_HIVADULT * NG * proj_years);
  deaths_paeds_art.attr("dim") = NumericVector::create(hDS, trans,  pIDX_HIVADULT, NG, proj_years);
  
  NumericVector grad_paeds(trans * hDS * pIDX_HIVADULT * NG * proj_years);
  grad_paeds.attr("dim") = NumericVector::create(hDS, trans,  pIDX_HIVADULT, NG, proj_years);
  
  NumericVector grad_paeds_art(hTS * hDS * pIDX_HIVADULT * NG * proj_years);
  grad_paeds_art.attr("dim") = NumericVector::create(hTS, hDS,  pIDX_HIVADULT, NG, proj_years);
  
  
  NumericVector init_art_paed(trans * hDS * pIDX_HIVADULT * NG * proj_years);
  init_art_paed.attr("dim") = NumericVector::create(hDS, trans,  pIDX_HIVADULT, NG, proj_years);
 
  if (hAG == hAG_FULL) {
    leapfrog_sim<double, NG, pAG, pIDX_FERT, pAG_FERT,
		 pIDX_HIVADULT, hAG_FULL, hDS, hDS_adol, trans, tx_time, hTS>
      (REAL(demp["basepop"]),
       REAL(demp["Sx"]),
       REAL(demp["netmigr_adj"]),
       REAL(demp["asfr"]),
       REAL(demp["births_sex_prop"]),
       REAL(projp["incidinput"]),
       REAL(projp["incrr_sex"]),
       REAL(projp["fert_rat"]),
       REAL(projp["incrr_age"]),
       REAL(projp["cd4_initdist_full"]),
       REAL(projp["cd4_prog_full"]),
       REAL(projp["cd4_mort_full"]),
       REAL(projp["art_mort_full"]),
       REAL(projp["artmx_timerr"]),
       REAL(projp["art15plus_num"]),
       REAL(projp["scalar_art"]),
       LOGICAL(projp["art15plus_isperc"]),
       INTEGER(projp["artcd4elig_idx"]),
       *INTEGER(projp["art_alloc_method"]),
       *REAL(projp["art_alloc_mxweight"]),
       *INTEGER(projp["scale_cd4_mort"]),
       REAL(projp["art_dropout"]),
       REAL(projp["age15hivpop"]),
       REAL(projp["paed_incid_input"]),
       REAL(projp["paed_cd4_dist"]),
       REAL(projp["paed_cd4_prog"]),
       REAL(projp["adol_cd4_prog"]),
       REAL(projp["paed_cd4_mort"]),
       REAL(projp["adol_cd4_mort"]),
       REAL(projp["paed_art_mort"]),
       REAL(projp["adol_art_mort"]),
       REAL(projp["ctx_val"]),
       REAL(projp["paed_cd4_transition"]),
       ctx_effect,
       REAL(projp["paed_art_val"]),
       LOGICAL(projp["artpaeds_isperc"]),
       REAL(projp["paed_art_elig_age"]),
       REAL(projp["paed_art_elig_cd4"]),
       proj_years,
       hiv_steps_per_year,
       *INTEGER(projp["t_ART_start"]) - 1, // 0-based indexing vs. R 1-based
       INTEGER(projp["hAG_SPAN_full"]),
       REAL(totpop1),
       REAL(hivpop1),
       REAL(hivnpop1),
       REAL(infections),
       REAL(hivstrat_adult),
       REAL(artstrat_adult),
       REAL(hivstrat_paeds),
       REAL(artstrat_paeds),
       REAL(artelig_paeds),
       REAL(births),
       REAL(hiv_births),
       REAL(natdeaths),
       REAL(natdeaths_hivpop),
       REAL(hivdeaths),
       REAL(aidsdeaths_noart),
       REAL(aidsdeaths_art),
       REAL(aidsdeaths_noart_paed),
       REAL(aidsdeaths_art_paed),
       REAL(artnum_paed),
       REAL(artinit),
       REAL(deaths_paeds),
       REAL(deaths_paeds_art),
       REAL(grad_paeds),
       REAL(grad_paeds_art),
       REAL(init_art_paed),
       REAL(coarse_totpop1));
  } else if (hAG == hAG_COARSE) {
    leapfrog_sim<double, NG, pAG, pIDX_FERT, pAG_FERT,
		 pIDX_HIVADULT, hAG_COARSE, hDS, hDS_adol, trans, tx_time,  hTS>
      (REAL(demp["basepop"]),
       REAL(demp["Sx"]),
       REAL(demp["netmigr_adj"]),
       REAL(demp["asfr"]),
       REAL(demp["births_sex_prop"]),
       REAL(projp["incidinput"]),
       REAL(projp["incrr_sex"]),
       REAL(projp["fert_rat"]),
       REAL(projp["incrr_age"]),
       REAL(projp["cd4_initdist_coarse"]),
       REAL(projp["cd4_prog_coarse"]),
       REAL(projp["cd4_mort_coarse"]),
       REAL(projp["art_mort_coarse"]),
       REAL(projp["artmx_timerr"]),
       REAL(projp["art15plus_num"]),
       REAL(projp["scalar_art"]),
       LOGICAL(projp["art15plus_isperc"]),
       INTEGER(projp["artcd4elig_idx"]),
       *INTEGER(projp["art_alloc_method"]),
       *REAL(projp["art_alloc_mxweight"]),
       *INTEGER(projp["scale_cd4_mort"]),
       REAL(projp["art_dropout"]), 
       REAL(projp["age15hivpop"]),
       REAL(projp["paed_incid_input"]),
       REAL(projp["paed_cd4_dist"]),
       REAL(projp["paed_cd4_prog"]),
       REAL(projp["adol_cd4_prog"]),
       REAL(projp["paed_cd4_mort"]),
       REAL(projp["adol_cd4_mort"]),
       REAL(projp["paed_art_mort"]),
       REAL(projp["adol_art_mort"]),
       REAL(projp["ctx_val"]),
       REAL(projp["paed_cd4_transition"]),
       ctx_effect,
       REAL(projp["paed_art_val"]),
       LOGICAL(projp["artpaeds_isperc"]),
       REAL(projp["paed_art_elig_age"]),
       REAL(projp["paed_art_elig_cd4"]),
       proj_years,
       hiv_steps_per_year,
       *INTEGER(projp["t_ART_start"]) - 1,  // 0-based indexing vs. R 1-based
       INTEGER(projp["hAG_SPAN_coarse"]),
       REAL(totpop1),
       REAL(hivpop1),
       REAL(hivnpop1),
       REAL(infections),
       REAL(hivstrat_adult),
       REAL(artstrat_adult),
       REAL(hivstrat_paeds),
       REAL(artstrat_paeds),
       REAL(artelig_paeds),
       REAL(births),
       REAL(hiv_births),
       REAL(natdeaths),
       REAL(natdeaths_hivpop),
       REAL(hivdeaths),
       REAL(aidsdeaths_noart),
       REAL(aidsdeaths_art),
       REAL(aidsdeaths_noart_paed),
       REAL(aidsdeaths_art_paed),
       REAL(artnum_paed),
       REAL(artinit),
       REAL(deaths_paeds),
       REAL(deaths_paeds_art),
       REAL(grad_paeds),
       REAL(grad_paeds_art),
       REAL(init_art_paed),
       REAL(coarse_totpop1));
  } else {
    Rf_error("Invalid HIV stratification age groups (hAG)");
  }

  
  List ret = ListBuilder()
    .add("totpop1", totpop1)
    .add("hivpop1", hivpop1)
    .add("hivnpop1", hivnpop1)
    .add("hivstrat_adult", hivstrat_adult)
    .add("artstrat_adult", artstrat_adult)
    .add("hivstrat_paeds", hivstrat_paeds)
    .add("artstrat_paeds", artstrat_paeds)
    .add("artelig_paeds", artelig_paeds)
    .add("infections", infections)
    .add("births", births)
    .add("hiv_births", hiv_births)
    .add("natdeaths", natdeaths)
    .add("natdeaths_hivpop", natdeaths_hivpop)
    .add("hivdeaths", hivdeaths)
    .add("aidsdeaths_noart", aidsdeaths_noart)
    .add("aidsdeaths_art", aidsdeaths_art)
    .add("aidsdeaths_noart_paed", aidsdeaths_noart_paed)
    .add("aidsdeaths_art_paed", aidsdeaths_art_paed)
    .add("artnum_paed", artnum_paed)
    .add("artinit", artinit)
    .add("deaths_paeds", deaths_paeds)
    .add("deaths_paeds_art", deaths_paeds_art)
    .add("grad_paeds", grad_paeds)
    .add("grad_paeds_art", grad_paeds_art)
    .add("init_art_paed", init_art_paed)
    .add("coarse_totpop1", coarse_totpop1);

  return ret;
}
