#ifndef LEAPFROG_H
#define LEAPFROG_H

#include <unsupported/Eigen/CXX11/Tensor>

using Eigen::TensorBase;
using Eigen::TensorMap;
using Eigen::Tensor;
using Eigen::TensorFixedSize;
using Eigen::Sizes;

//' The output memory is passed as an argument rather than constructed
//' by the function.
//'
//' Template parameters
//' @param Type data type variables involved in calculations (typically 'double'
//'   but required to be a templated parameter for autodiff integration).
//' @param NG number of genders (= 2)
//' @param pAG number of population age groups (e.g. 81 for ages 0 to 80+)
//' @param pIDX_FERT first index eligible for fertility
//' @param pAG_FERT number of ages eligible for fertility
//' @param hAG number of age categories in the adult population
//' @param hDS number of CD4 categories
//' @param hTS length on ART (less than 6 months, 6-12 months, 12+ months)
//' @param hTM transmission category
//' @param hPS number of PMTCT types (7)

//'
//' @details
//' State space dimensions are specified as template parameters so that
//' these are specified at compile time, allowing stack allocation of
//' working arrays.
//'
//' Notes:
//' * Consider adding half the migrants at the start and half at the end?
//' * Unsure whether to have this function accept raw pointers or TensorMap
//'   - If passing Tensor objects, it would be ideal to pass TensorBase, but
//'     this is not possible because dimensions are not known.

template <typename Type, int NG, int pAG, int pIDX_FERT, int pAG_FERT,
  int pIDX_HIVADULT,
  int hAG, int hDS, int hDS_adol,int hTM, int hTS, int hPS, int hBF>
  void leapfrog_sim(const Type *p_basepop,
                    const Type *p_sx,
                    const Type *p_netmigr,
                    const Type *p_asfr,
                    const Type *p_births_sex_prop,
                    //
                    // adult incidence
                    const Type *p_incidinput,
                    const Type *p_incrr_sex,
                    const Type *p_tfr,
                    const Type *p_incrr_age,
                    //
                    // adult natural history
                    const Type *p_cd4_initdist,
                    const Type *p_cd4_prog,
                    const Type *p_cd4_mort,
                    const Type *p_art_mort,
                    const Type *p_artmx_timerr,
                    //
                    // adult ART
                    const Type *p_art15plus_num,
                    const int *p_art15plus_isperc,
                    const int *artcd4elig_idx,
                    const int art_alloc_method,
                    const Type art_alloc_mxweight,
                    const int scale_cd4_mort,
                    const Type *p_art_dropout,
                    const Type *p_age15hivpop,

                    //paed inputs
                    const Type *p_paed_incid_input,
                    const Type *p_paed_cd4_dist,
                    const Type *p_paed_cd4_prog,
                    const Type *p_adol_cd4_prog,
                    const Type *p_paed_cd4_mort,
                    const Type *p_adol_cd4_mort,
                    const Type *p_paed_art_mort,
                    const Type *p_adol_art_mort,
                    const Type *p_ctx_val,
                    const Type *p_init_art_dist,
                    const Type *p_mort_art_rr,
                    const Type *p_paed_cd4_transition,
                    const Type *p_adult_cd4_dist,
                    const double ctx_effect,
                    const Type *p_paed_art_val,
                    const int *p_artpaeds_isperc,
                    const int *p_pmtct_input_isperc,
                    const Type *p_paed_art_elig_age,
                    const Type *p_paed_art_elig_cd4,
                    const Type *p_paed_art_ltfu,

                    //mtct trans
                    const Type *p_mtct_trans,
                    const Type *p_fert_mult_by_age,
                    const Type *p_fert_mult_offart,
                    const Type *p_fert_mult_onart,
                    const Type *p_pmtct_mtct,
                    const Type *p_pmtct,
                    const Type *p_bf_duration,
                    const Type *p_pmtct_dropout,
                    //
                    //settings
                    const int sim_years,
                    const int hiv_steps_per_year,
                    const int t_ART_start,
                    const int *hAG_SPAN,
                    //
                    //outputs
                    Type *p_totpop1,
                    Type *p_hivpop1,
                    Type *p_hivnpop1,
                    Type *p_infections,
                    Type *p_hivstrat_adult,
                    Type *p_artstrat_adult,
                    Type *p_hivstrat_paeds,
                    Type *p_artstrat_paeds,
                    Type *p_artelig_paeds,
		                Type *p_births,
		                Type *p_hiv_births,
                    Type *p_natdeaths,
                    Type *p_natdeaths_hivpop,
                    Type *p_hivdeaths,
                    Type *p_aidsdeaths_noart,
                    Type *p_aidsdeaths_art,
                    Type *p_aidsdeaths_noart_paed,
                    Type *p_aidsdeaths_art_paed,
                    Type *p_artnum_paed,
                    Type *p_artinit,
                    Type *p_deaths_paeds,
                    Type *p_deaths_paeds_art,
                    Type *p_grad_paeds,
                    Type *p_grad_paeds_art,
                    Type *p_init_art_paed,
                    Type *p_coarse_totpop1,
                    Type *p_tracking) {



  // macros
  // TODO: unsure if these should be defined here or elsewhere. Maybe in a future class.
  typedef Eigen::TensorMap<Eigen::Tensor<const int, 1>> TensorMapX1cI;
  typedef Eigen::TensorMap<Eigen::Tensor<const int, 2>> TensorMapX2cI;

  typedef Eigen::TensorMap<Eigen::Tensor<const Type, 1>> TensorMapX1cT;
  typedef Eigen::TensorMap<Eigen::Tensor<const Type, 2>> TensorMapX2cT;
  typedef Eigen::TensorMap<Eigen::Tensor<const Type, 3>> TensorMapX3cT;
  typedef Eigen::TensorMap<Eigen::Tensor<const Type, 4>> TensorMapX4cT;

  typedef Eigen::TensorMap<Eigen::Tensor<Type, 1>> TensorMapX1T;
  typedef Eigen::TensorMap<Eigen::Tensor<Type, 2>> TensorMapX2T;
  typedef Eigen::TensorMap<Eigen::Tensor<Type, 3>> TensorMapX3T;
  typedef Eigen::TensorMap<Eigen::Tensor<Type, 4>> TensorMapX4T;
  typedef Eigen::TensorMap<Eigen::Tensor<Type, 5>> TensorMapX5T;

  const int MALE = 0;
  const int FEMALE = 1;

  const int ART0MOS = 0;

  const int pIDX_INCIDPOP = pIDX_HIVADULT;
  const int pAG_INCIDPOP = 35; // !!! HARD CODED 15-49 for now

  const int hAG_15PLUS = hAG;
  const int hIDX_15PLUS = 0;
  const Type h_art_stage_dur[hTS-1] = {0.5, 0.5};


  // // inputs

  // state space dimensions


  double dt = 1.0 / hiv_steps_per_year;

  // demography
  const TensorMapX3cT basepop(p_basepop, pAG, NG, sim_years);
  const TensorMapX3cT sx(p_sx, pAG+1, NG, sim_years);
  const TensorMapX3cT netmigr(p_netmigr, pAG, NG, sim_years);
  const TensorMapX2cT asfr(p_asfr, pAG_FERT, sim_years);
  const TensorMapX2cT births_sex_prop(p_births_sex_prop, NG, sim_years);

  // adult HIV incidence
  const TensorMapX1cT incidinput(p_incidinput, sim_years);
  const TensorMapX2cT incrr_sex(p_incrr_sex, NG, sim_years);
  const TensorMapX1cT tfr(p_tfr, sim_years);
  const TensorMapX3cT incrr_age(p_incrr_age, pAG - pIDX_HIVADULT, NG, sim_years);

  // adult HIV natural history
  const TensorMapX3cT cd4_initdist(p_cd4_initdist, hDS, hAG, NG);
  const TensorMapX3cT cd4_prog(p_cd4_prog, hDS-1, hAG, NG);
  const TensorMapX3cT cd4_mort(p_cd4_mort, hDS, hAG, NG);
  const TensorMapX4cT art_mort(p_art_mort, hTS, hDS, hAG, NG);
  const TensorMapX2cT artmx_timerr(p_artmx_timerr, hTS, sim_years);

  // adult ART
  const TensorMapX2cT art15plus_num(p_art15plus_num, NG, sim_years);
  const TensorMapX2cI art15plus_isperc(p_art15plus_isperc, NG, sim_years);
  const TensorMapX1cT art_dropout(p_art_dropout, sim_years);
  const TensorMapX4cT age15hivpop(p_age15hivpop, 4, hDS, NG, sim_years);

  // paed
  const TensorMapX1cT paed_incid_input(p_paed_incid_input, sim_years);
  const TensorMapX1cT paed_cd4_dist(p_paed_cd4_dist, hDS);
  const TensorMapX1cT paed_cd4_prog(p_paed_cd4_prog, hDS);
  const TensorMapX1cT adol_cd4_prog(p_adol_cd4_prog, hDS_adol);
  const TensorMapX3cT paed_cd4_mort(p_paed_cd4_mort, hDS, hTM, 5);
  const TensorMapX3cT adol_cd4_mort(p_adol_cd4_mort, hDS_adol, hTM, 10);
  const TensorMapX1cT ctx_val(p_ctx_val, sim_years);
  const TensorMapX1cT paed_art_val(p_paed_art_val, sim_years);
  const TensorMapX1cT paed_art_elig_age(p_paed_art_elig_age, sim_years);
  const TensorMapX1cT paed_art_ltfu(p_paed_art_ltfu, sim_years);
  const TensorMapX2cT init_art_dist(p_init_art_dist, pIDX_HIVADULT, sim_years);
  const TensorMapX3cT mort_art_rr(p_mort_art_rr, hTS, pIDX_HIVADULT, sim_years);



  //5 age categories
  const TensorMapX2cT paed_art_elig_cd4(p_paed_art_elig_cd4, pIDX_FERT, sim_years);
  //count by pct
  const TensorMapX2cT paed_cd4_transition(p_paed_cd4_transition, hDS_adol, hDS);
  const TensorMapX2cT adult_cd4_dist(p_adult_cd4_dist, hDS, hDS_adol);

  const TensorMapX3cT paed_art_mort(p_paed_art_mort, hDS, hTS, 5);
  const TensorMapX3cT adol_art_mort(p_adol_art_mort, hDS_adol, hTS, 10);
  const TensorMapX1cI artpaeds_isperc(p_artpaeds_isperc, sim_years);
  const TensorMapX1cI pmtct_input_isperc(p_pmtct_input_isperc, sim_years);


  //MTCT
  const TensorMapX1cT mtct_trans(p_mtct_trans, hDS);
  const TensorMapX1cT fert_mult_by_age(p_fert_mult_by_age, 35);
  const TensorMapX1cT fert_mult_offart(p_fert_mult_offart, hDS);
  const TensorMapX1cT fert_mult_onart(p_fert_mult_onart, 35);
  //2 is perinatal and bf
  //3 is ART <4 weeks, ART >4 weeks, ART before
  //2 is the same as above, but 5 is no trt, A, B, nev, WHO 2006
  const TensorMapX3cT pmtct_mtct(p_pmtct_mtct, hDS, 8, 2);
  //hPS is number of pmtct types (7) and two is whether it is a number or percent
  const TensorMapX3cT pmtct(p_pmtct, hPS, sim_years, 2);
  //18 is the number of 2 month age groups up to 3 years
  const TensorMapX3cT bf_duration(p_bf_duration, 18, sim_years, 2);
  //6 is for percent already on ART retained at delivery
  //percent starting ART retained at delivery
  //postnatal option A monthly dropout
  //postnatal option B monthly dropout
  //art 0-12 months breastfeeding
  //art 12+months bf
  const TensorMapX2cT pmtct_dropout(p_pmtct_dropout, 6, sim_years);

  // outputs
  TensorMapX3T totpop1(p_totpop1, pAG, NG, sim_years);
  TensorMapX3T hivpop1(p_hivpop1, pAG, NG, sim_years);
  TensorMapX3T hivnpop1(p_hivnpop1, pAG, NG, sim_years);
  TensorMapX3T infections(p_infections, pAG, NG, sim_years);
  TensorMapX1T births(p_births, sim_years);
  TensorMapX1T hiv_births(p_hiv_births,  sim_years);
  //TensorMapX2T hiv_births(p_hiv_births, hDS, sim_years);

  TensorMapX3T natdeaths(p_natdeaths, pAG, NG, sim_years);
  TensorMapX3T natdeaths_hivpop(p_natdeaths_hivpop, pAG, NG, sim_years);
  TensorMapX3T hivdeaths(p_hivdeaths, pAG, NG, sim_years);
  TensorMapX4T aidsdeaths_noart(p_aidsdeaths_noart, hDS, hAG,  NG, sim_years);
  TensorMapX5T aidsdeaths_art(p_aidsdeaths_art, hTS, hDS, hAG, NG, sim_years);
  TensorMapX4T artinit(p_artinit, hDS, hAG, NG, sim_years);

  TensorMapX4T hivstrat_adult(p_hivstrat_adult, hDS, hAG, NG, sim_years);
  TensorMapX5T artstrat_adult(p_artstrat_adult, hTS, hDS, hAG, NG, sim_years);
  //adding in transmission category
  TensorMapX5T hivstrat_paeds(p_hivstrat_paeds, hDS, hTM, pIDX_HIVADULT, NG, sim_years);
  TensorMapX5T artstrat_paeds(p_artstrat_paeds, hTS, hDS, pIDX_HIVADULT, NG, sim_years);
  TensorMapX5T artelig_paeds(p_artelig_paeds, hDS, hTM, pIDX_HIVADULT, NG, sim_years);
  TensorMapX5T aidsdeaths_noart_paed(p_aidsdeaths_noart_paed, hDS, hTM, pIDX_HIVADULT, NG, sim_years);
  TensorMapX5T aidsdeaths_art_paed(p_aidsdeaths_art_paed, hTS, hDS, pIDX_HIVADULT, NG, sim_years);
  TensorMapX5T deaths_paeds(p_deaths_paeds, hDS, hTM, pIDX_HIVADULT, NG, sim_years);
  TensorMapX5T deaths_paeds_art(p_deaths_paeds_art, hDS, hTM, pIDX_HIVADULT, NG, sim_years);
  TensorMapX5T grad_paeds(p_grad_paeds, hDS, hTM, pIDX_HIVADULT, NG, sim_years);
  TensorMapX5T grad_paeds_art(p_grad_paeds_art, hTS, hDS,  pIDX_HIVADULT, NG, sim_years);
  TensorMapX1T artnum_paed(p_artnum_paed, sim_years);
  TensorMapX5T init_art_paed(p_init_art_paed, hDS, hTM, pIDX_HIVADULT, NG, sim_years);


  TensorMapX3T coarse_totpop1(p_coarse_totpop1, hAG, NG, sim_years);
  TensorMapX3T tracking(p_tracking, 18, 35, sim_years);





  // initialise population

  for(int g = 0; g < NG; g++) {
    for(int a = 0; a < pAG; a++) {
      totpop1(a, g, 0) = basepop(a, g, 0);
    }
  }

  hivpop1.setZero();
  hivstrat_adult.setZero();
  artstrat_adult.setZero();
  hivstrat_paeds.setZero();
  artstrat_paeds.setZero();
  artelig_paeds.setZero();
  aidsdeaths_noart_paed.setZero();
  aidsdeaths_art_paed.setZero();

  deaths_paeds.setZero();
  deaths_paeds_art.setZero();
  grad_paeds.setZero();
  grad_paeds_art.setZero();
  init_art_paed.setZero();
  hiv_births.setZero();
  tracking.setZero();

  int everARTelig_idx = hDS;


  /////ISSUE: doesn't like this sim_years call for some reason?
  //TensorFixedSize<Type, Sizes<sim_years>> needPMTCT;
  TensorFixedSize<Type, Sizes<61>> needPMTCT;
  TensorFixedSize<Type, Sizes<61>> HIVPregwomen;
  TensorFixedSize<Type, Sizes<61>> OnPMTCT;




  ////////////////////////////////////
  ////  do population projection  ////
  ////////////////////////////////////

  for(int t = 1; t < sim_years; t++){




    TensorFixedSize<Type, Sizes< hDS, NG>> age15_hivpop;
    age15_hivpop.setZero();
    for(int g = 0; g < NG; g++){
      for(int hm = 0; hm < hDS; hm++){
        for(int cat = 0; cat < 4; cat++){
          age15_hivpop(hm, g) += hivstrat_paeds(hm, cat, 14, g, t-1) ;

        }
      }
    }



    TensorFixedSize<Type, Sizes< hTS, hDS, NG>> age15_artpop;
    age15_artpop.setZero();
    for(int g = 0; g < NG; g++){
      for(int hm = 0; hm < hDS; hm++){
        for(int hu = 0; hu < hTS; hu ++){
          //if(art15plus_num(g,t) > 0){
            age15_artpop(hu, hm, g) = artstrat_paeds(hu, hm, 14, g, t-1) ;
        //  }
        }
      }
    }


    if(t == 19){
      double mag;
      mag = 0.0;
      for(int g = 0; g < NG; g++){
        for(int hm = 0; hm < hDS; hm++){
            mag += age15_hivpop(hm, g)  ;
        }
      }
   //   std::cout << mag;
    }

    TensorFixedSize<Type, Sizes<pAG, NG>> migrate_ag;
    // ageing and non-HIV mortality
    for(int g = 0; g < NG; g++){

      for(int a = 1; a < pAG; a++) {
        natdeaths(a, g, t) = totpop1(a-1, g, t-1) * (1.0 - sx(a, g, t));
        totpop1(a, g, t) = totpop1(a-1, g, t-1) - natdeaths(a, g, t);
      }

      // open age group
      Type natdeaths_open_age = totpop1(pAG-1, g, t-1) * (1.0 - sx(pAG, g, t));
      natdeaths(pAG-1, g, t) += natdeaths_open_age;
      totpop1(pAG-1, g, t) += totpop1(pAG-1, g, t-1) - natdeaths_open_age;

      // net migration
      for(int a = 1; a < pAG - 1; a++) {
	// Number of net migrants adjusted for survivorship to end of period (qx / 2)
        migrate_ag(a, g) = netmigr(a, g, t) * (1.0 + sx(a, g, t)) * 0.5 / totpop1(a, g, t);
        totpop1(a, g, t) *= 1.0 + migrate_ag(a, g);
      }

      // For open age group, netmigrant survivor adjustment based on weighted
      // sx for age 79 and age 80+.
      // * Numerator: totpop1(a, g, t-1) * (1.0 + sx(a+1, g, t)) + totpop1(a-1, g, t-1) * (1.0 + sx(a, g, t))
      // * Denominator: totpop1(a, g, t-1) + totpop1(a-1, g, t-1)
      // Re-expressed current population and deaths to open age group (already calculated):
      int a = pAG - 1;
      Type sx_netmig = (totpop1(a, g,t) + 0.5 * natdeaths(pAG-1, g, t)) / (totpop1(a, g,t) + natdeaths(pAG-1, g, t));
      migrate_ag(a, g) = sx_netmig * netmigr(a, g, t) / totpop1(a, g,t);
      totpop1(a, g, t) *= 1.0 + migrate_ag(a, g);
    }



    // // demographic projection of the adult HIV population

    // age population and calculate non-HIV deaths to HIV population
    for(int g = 0; g < NG; g++){
      for(int a = 1; a < pAG; a++) {
        natdeaths_hivpop(a, g, t) = hivpop1(a-1, g, t-1) * (1.0 - sx(a, g, t));
        hivpop1(a, g, t) = hivpop1(a-1, g, t-1);
      }

      // open age group
      natdeaths_hivpop(pAG-1, g, t) += hivpop1(pAG-1, g, t-1) * (1.0 - sx(pAG, g, t));
      hivpop1(pAG-1, g, t) += hivpop1(pAG-1, g, t-1);
    }




    // age coarse stratified HIV population
    TensorFixedSize<Type, Sizes<hAG, NG>> hiv_ag_prob;
    hiv_ag_prob.setZero();
    TensorFixedSize<Type, Sizes<hAG, NG>> tot_ag_prob;
    tot_ag_prob.setZero();
    TensorFixedSize<Type, Sizes<pIDX_FERT, NG>> hivpaeds_ag_prob;
    hivpaeds_ag_prob.setZero();


    for(int g = 0; g < NG; g++){
      int a = pIDX_HIVADULT ;
      for(int ha = 0; ha < (hAG-1); ha++){
        for(int i = 0; i < hAG_SPAN[ha]; i++){
          hiv_ag_prob(ha, g) += hivpop1(a, g, t-1);
          tot_ag_prob(ha, g) += totpop1(a, g, t-1);
          //save tot_ag_prob out here as the first index of coarse_totpop1
          a++;
        }
        hiv_ag_prob(ha, g) = (hiv_ag_prob(ha, g) > 0) ? hivpop1(a-1, g, t-1) / hiv_ag_prob(ha, g) : 0.0;
        tot_ag_prob(ha, g) = (tot_ag_prob(ha, g) > 0) ? totpop1(a-1, g, t-1) / tot_ag_prob(ha, g) : 0.0;

      }
      // Note: loop stops at hAG-1; no one ages out of the open-ended age group
    }

    for(int g = 0; g < NG; g++) {
      for(int ha = 1; ha < hAG; ha++) {
        coarse_totpop1(ha, g, t) = (1.0 - tot_ag_prob(ha, g)) * coarse_totpop1(ha, g, t-1);  // age-out
        coarse_totpop1(ha, g, t) += tot_ag_prob(ha-1, g) * coarse_totpop1(ha-1, g, t-1);   // age-in
        for(int hm = 0; hm < hDS; hm++) {
          hivstrat_adult(hm, ha, g, t) = (1.0 - hiv_ag_prob(ha, g)) * hivstrat_adult(hm, ha, g, t-1);  // age-out
          hivstrat_adult(hm, ha, g, t) += hiv_ag_prob(ha-1, g) * hivstrat_adult(hm, ha-1, g, t-1);   // age-in
          //MAGGIE: This needs to be changed, ART issue
          if(t > t_ART_start or paed_art_val(t-1) > 0)
            for(int hu = 0; hu < hTS; hu++) {
              artstrat_adult(hu, hm, ha, g, t) = (1.0 - hiv_ag_prob(ha, g)) * artstrat_adult(hu, hm, ha, g, t-1);
              artstrat_adult(hu, hm, ha, g, t) += hiv_ag_prob(ha-1, g) * artstrat_adult(hu, hm, ha-1, g, t-1);
            }
        }
      }
    }


    // !!!TODO: add HIV+ 15 year old entrants
    for (int g = 0; g < NG; g++) {
      for (int hm = 0; hm < hDS; hm++) {
        for(int hm_adol = 0; hm_adol < hDS_adol; hm_adol++){
          hivstrat_adult(hm, 0, g, t) += age15_hivpop(hm_adol, g) * adult_cd4_dist(hm, hm_adol) ;//* sx(pIDX_HIVADULT,g,t) ;
          for(int hu = 0; hu < hTS; hu++) {
          artstrat_adult(hu, hm, 0, g, t) += age15_artpop(hu, hm_adol, g) * adult_cd4_dist(hm, hm_adol);// * sx(pIDX_HIVADULT,g,t) ;
          }
        }
      }
    }








    for(int g = 0; g < NG; g++){
      for(int ha = 1; ha < 5; ha++) {
        for(int hm = 0; hm < hDS; hm++){
          for(int cat = 0 ; cat < hTM; cat++){
           hivstrat_paeds(hm, cat, ha, g, t) += hivstrat_paeds(hm, cat, ha-1, g, t-1) * sx(ha, g, t);

          }
            for(int dur = 0; dur < hTS; dur++){
              artstrat_paeds(dur, hm, ha, g, t) += artstrat_paeds(dur, hm, ha-1, g, t-1) * sx(ha, g, t);
            }

        }
      }
    }



    for(int g = 0; g < NG; g++){
        for(int hm = 0; hm < hDS; hm++){
          for(int hm_alt = 0; hm_alt < hDS_adol; hm_alt++){
          for(int cat = 0 ; cat < hTM; cat++){
             hivstrat_paeds(hm_alt, cat, 5, g, t) +=  hivstrat_paeds(hm, cat, 4, g, t-1) * sx(5, g, t) * paed_cd4_transition(hm_alt, hm);
        }
          for(int dur = 0; dur < hTS; dur++){
            artstrat_paeds(dur, hm_alt, 5, g, t) += artstrat_paeds(dur, hm, 4, g, t-1) * sx(5, g, t) * paed_cd4_transition(hm_alt, hm);
          }
      }
    }
  }




    for(int g = 0; g < NG; g++){
      for(int ha = 6; ha < pIDX_HIVADULT; ha++) {
        for(int hm = 0; hm < hDS_adol; hm++){
          for(int cat = 0 ; cat < hTM; cat++){
           hivstrat_paeds(hm, cat, ha, g, t) += hivstrat_paeds(hm, cat, ha-1, g, t-1) * sx(ha, g, t);

          }
          for(int dur = 0; dur < hTS; dur++){
            artstrat_paeds(dur, hm, ha, g, t) += artstrat_paeds(dur, hm, ha-1, g, t-1) * sx(ha, g, t);
          }

        }
      }
    }



    TensorFixedSize<Type, Sizes<hAG, NG>> hivpop_ha;
    hivpop_ha.setZero();
    for(int g = 0; g < NG; g++) {
      int a = pIDX_HIVADULT;
      for(int ha = 0; ha < hAG; ha++){
        for(int i = 0; i < hAG_SPAN[ha]; i++){
          hivpop_ha(ha, g) += hivpop1(a, g, t);
          a++;
        }
      }
    }


    // remove non-HIV deaths and net migration from hivpop1
    TensorFixedSize<Type, Sizes<pAG, NG>> netmig_ag;
    for(int g = 0; g < NG; g++) {
      for(int a = 1; a < pAG; a++) {
        hivpop1(a, g, t) -= natdeaths_hivpop(a, g, t);
        netmig_ag(a, g) = hivpop1(a, g, t) * migrate_ag(a, g);
        hivpop1(a, g, t) += netmig_ag(a, g);
      }
    }


    // remove non-HIV deaths and net migration from adult stratified population
    for(int g = 0; g < NG; g++){
      int a = pIDX_HIVADULT;
      for(int ha = 0; ha < hAG; ha++){
        Type deathsmig_ha = 0;
        for(int i = 0; i < hAG_SPAN[ha]; i++){
          deathsmig_ha -= natdeaths_hivpop(a, g, t);
          deathsmig_ha += netmig_ag(a, g);
          a++;
        }

        //implement this for children
        Type deathmigrate_ha = hivpop_ha(ha, g) > 0 ? deathsmig_ha / hivpop_ha(ha, g) : 0.0;
        for(int hm = 0; hm < hDS; hm++){
          hivstrat_adult(hm, ha, g, t) *= 1.0 + deathmigrate_ha;
          //MAGGIE: This needs to be changed, ART issue
          if(t > t_ART_start or paed_art_val(t-1) > 0) {
            for(int hu = 0; hu < hTS; hu++) {
              artstrat_adult(hu, hm, ha, g, t) *= 1.0 + deathmigrate_ha;
            }
          }
        } // loop over hm
      } // loop over ha
    } // loop over g



    // fertility


    births(t) = 0.0;
   for(int af = 0; af < pAG_FERT; af++) {
     births(t) += (totpop1(pIDX_FERT + af, FEMALE, t-1) + totpop1(pIDX_FERT + af, FEMALE, t)) * 0.5 * asfr(af, t);

    }


    // add births
    for(int g = 0; g < NG; g++) {
      Type births_sex = births(t) * births_sex_prop(g, t);
      natdeaths(0, g, t) = births_sex * (1.0 - sx(0, g, t));
      totpop1(0, g, t) =  births_sex * sx(0, g, t);

      // Assume 2/3 survival rate since mortality in first six months higher than
      // second 6 months (Spectrum manual, section 6.2.7.4)
      Type migrate_a0 = netmigr(0, g, t) * (1.0 + 2.0 * sx(0, g, t)) / 3.0 / totpop1(0, g, t);
      totpop1(0, g, t) *= 1.0 + migrate_a0;
    }


    //////////////////////////////////////
    ////  Adult HIV model simulation  ////
    //////////////////////////////////////

    int cd4elig_idx = artcd4elig_idx[t] - 1; // -1 for 0-based indexing vs. 1-based in R
    everARTelig_idx = cd4elig_idx < everARTelig_idx ? cd4elig_idx : everARTelig_idx;
    // !! SPECIAL POPULATION ART ELIGIBILITY NOT YET IMPLEMENTED;
    //    will replace previous line when implemented
    // int anyelig_idx = (specpop_percelig[t] > 0 | pw_artelig[t] > 0) ? 0 : cd4elig_idx;
    // everARTelig_idx = anyelig_idx < everARTelig_idx ? anyelig_idx : everARTelig_idx;
    int anyelig_idx = cd4elig_idx;

    TensorFixedSize<Type, Sizes<pAG, NG>> infections_ts;

    const int EPP_DIRECTINCID_HTS = 0;
    const int EPP_DIRECTINCID_ANN = 1;
    const int eppmod = EPP_DIRECTINCID_HTS;

    Type Xhivn[NG], incrate_g[NG];
    if(eppmod == EPP_DIRECTINCID_ANN ||
       eppmod == EPP_DIRECTINCID_HTS ) {

      // Calculate incidence rate by sex
      //
      // Note: In Spectrum, incidence rate by age is calculated per time-step using the
      // **current year** HIV population, instead of the previous year HIV population.
      // Rob Glaubius, 5 August 2022: https://github.com/mrc-ide/leapfrog/issues/18

      for(int g = 0; g < NG; g++){
        Xhivn[g] = 0.0;
        for(int a = pIDX_INCIDPOP; a < pIDX_INCIDPOP + pAG_INCIDPOP; a++){
          Xhivn[g] += totpop1(a, g, t-1) - hivpop1(a, g, t-1);
        }
      }


      Type incrate_i = incidinput(t);
      incrate_g[MALE] = incrate_i * (Xhivn[MALE]+Xhivn[FEMALE]) / (Xhivn[MALE] + incrr_sex(t)*Xhivn[FEMALE]);
      incrate_g[FEMALE] = incrate_i * incrr_sex(t)*(Xhivn[MALE]+Xhivn[FEMALE]) / (Xhivn[MALE] + incrr_sex(t)*Xhivn[FEMALE]);
    }


    for(int hts = 0; hts < hiv_steps_per_year; hts++) {

      TensorFixedSize<Type, Sizes<hAG, NG>> hivdeaths_ha;
      hivdeaths_ha.setZero();

      // untreated population

      // disease progression and mortality
      TensorFixedSize<Type, Sizes<hDS, hAG, NG>> grad;
      grad.setZero();

      for(int g = 0; g < NG; g++) {
        for(int ha = 0; ha < hAG; ha++) {
          for(int hm = 0; hm < hDS; hm++) {

            Type cd4mx_scale = 1.0;
            // !!NOTE: Mortality scaling not yet implemented
            if (scale_cd4_mort &
		(t >= t_ART_start) &
                (hm >= everARTelig_idx) &
                (hivstrat_adult(hm, ha, g, t) > 0.0)) {
              Type artpop_hahm = 0.0;
              for(int hu = 0; hu < hTS; hu++) {
                artpop_hahm += artstrat_adult(hu, hm, ha, g, t);
              }
              cd4mx_scale = hivstrat_adult(hm, ha, g, t) / (hivstrat_adult(hm, ha, g, t) + artpop_hahm);
            }

            Type deaths = cd4mx_scale * cd4_mort(hm, ha, g) * hivstrat_adult(hm, ha, g, t);
            hivdeaths_ha(ha, g) += dt*deaths;
            aidsdeaths_noart(hm, ha, g, t) += dt*deaths;
            grad(hm, ha, g) = -deaths;
          }

          for(int hm = 1; hm < hDS; hm++) {
            grad(hm-1, ha, g) -= cd4_prog(hm-1, ha, g) * hivstrat_adult(hm-1, ha, g, t);
            grad(hm, ha, g) += cd4_prog(hm-1, ha, g) * hivstrat_adult(hm-1, ha, g, t);
          }
        }
      }

      // calculate new infections by age
      if (eppmod == EPP_DIRECTINCID_HTS) {

	for(int g = 0; g < NG; g++) {

	  TensorFixedSize<Type, Sizes<pAG, NG>> hivn_a;
	  for(int a = pIDX_INCIDPOP; a < pAG; a++) {
	    hivn_a(a) = totpop1(a, g, t) - hivpop1(a, g, t);
	  }

	  Type Xhivn_incagerr = 0.0;
	  for(int a = pIDX_INCIDPOP; a < pIDX_INCIDPOP + pAG_INCIDPOP; a++){
	    Xhivn_incagerr += incrr_age(a - pIDX_INCIDPOP, g, t) * hivn_a(a);
	  }

	  for(int a = pIDX_HIVADULT; a < pAG; a++) {
	    infections_ts(a, g) = hivn_a(a) * incrate_g[g] * incrr_age(a - pIDX_INCIDPOP, g, t) * Xhivn[g] / Xhivn_incagerr;
	  }
	}
      }

      // add new infections to HIV population
      for(int g = 0; g < NG; g++){
        int a = pIDX_HIVADULT;
        for(int ha = 0; ha < hAG; ha++){
          Type infections_a, infections_ha = 0.0;
          for(int i = 0; i < hAG_SPAN[ha]; i++){
            infections_ha += infections_a = infections_ts(a, g);
            infections(a, g, t) += dt*infections_a;
            hivpop1(a, g, t) += dt*infections_a;
            a++;
          }

          // add infections to grad hivpop
          for(int hm = 0; hm < hDS; hm++) {
            grad(hm, ha, g) += infections_ha * cd4_initdist(hm, ha, g);
          }
        }
      }

      // // ART progression, mortality, and initiation
      //MAGGIE: This needs to be changed, ART issue, JUST NOT FOR INITIATION

      if(t >= t_ART_start){

	TensorFixedSize<Type, Sizes<hTS, hDS, hAG, NG>> gradART;

        // progression and mortality
        for(int g = 0; g < NG; g++)
          for(int ha = 0; ha < hAG; ha++)
            for(int hm = everARTelig_idx; hm < hDS; hm++){

              for(int hu = 0; hu < hTS; hu++){
		Type deaths = art_mort(hu, hm, ha, g) * artmx_timerr(hu, t) * artstrat_adult(hu, hm, ha, g, t);
		hivdeaths_ha(ha, g) += dt * deaths;
                aidsdeaths_art(hu, hm, ha, g, t) += dt * deaths;
                gradART(hu, hm, ha, g) = -deaths;
              }

              for(int hu = 0; hu < (hTS - 1); hu++) {
                gradART(hu, hm, ha, g) += -artstrat_adult(hu, hm, ha, g, t) / h_art_stage_dur[hu];
                gradART(hu+1, hm, ha, g) += artstrat_adult(hu, hm, ha, g, t) / h_art_stage_dur[hu];
              }

              // ART dropout
              if (art_dropout(t) > 0) {
                for (int hu = 0; hu < hTS; hu++) {
                  grad(hm, ha, g) += art_dropout(t) * artstrat_adult(hu, hm, ha, g, t);
                  gradART(hu, hm, ha, g) -= art_dropout(t) * artstrat_adult(hu, hm, ha, g, t);
                }
              }

            }


        // ART initiation
        for(int g = 0; g < NG; g++){

          TensorFixedSize<Type, Sizes<hDS, hAG_15PLUS>> artelig_hahm;

          Type Xart_15plus = 0.0, Xartelig_15plus = 0.0, expect_mort_artelig15plus = 0.0;
          for(int ha = hIDX_15PLUS; ha < hAG; ha++) {
            for(int hm = everARTelig_idx; hm < hDS; hm++) {
              if(hm >= anyelig_idx){
                // Type prop_elig = (hm >= cd4elig_idx) ? 1.0 : specpop_percelig[t];
                Type prop_elig = 1.0;  // !!! TODO: implement special population ART eligibility
                Xartelig_15plus += artelig_hahm(hm, ha-hIDX_15PLUS) = prop_elig * hivstrat_adult(hm, ha, g, t);
                expect_mort_artelig15plus += cd4_mort(hm, ha, g) * artelig_hahm(hm, ha-hIDX_15PLUS);
              }
              for(int hu = 0; hu < hTS; hu++)
                Xart_15plus += artstrat_adult(hu, hm, ha, g, t) + dt * gradART(hu, hm, ha, g);
            }

            // // if pw_artelig, add pregnant women to artelig_hahm population
            // if(g == FEMALE & pw_artelig[t] > 0 & ha < hAG_FERT){
            //   Type frr_pop_ha = 0;
            //   for(int a =  hAG_START[ha]; a < hAG_START[ha]+hAG_SPAN[ha]; a++)
            //     frr_pop_ha += pop_t(a, g, HIVN, t); // add HIV- population
            //   for(int hm = 0; hm < hDS; hm++){
            //     frr_pop_ha += frr_cd4(hm, ha-hIDX_FERT, t) * hivstrat_adult(hm, ha, g, t);
            //     for(int hu = 0; hu < hTS; hu++)
            //       frr_pop_ha += frr_art(hu, hm, ha-hIDX_FERT, t) * artstrat_adult(hu, hm, ha, g, t);
            //   }
            //   for(int hm = anyelig_idx; hm < cd4elig_idx; hm++){
            //     Type pw_elig_hahm = births_by_ha[ha-hIDX_FERT] * frr_cd4(hm, ha-hIDX_FERT, t) * hivstrat_adult(hm, ha, g, t) / frr_pop_ha;
            //     artelig_hahm(hm, ha-hIDX_15PLUS) += pw_elig_hahm;
            //     Xartelig_15plus += pw_elig_hahm;
            //     expect_mort_artelig15plus += cd4_mort(hm, ha, g) * pw_elig_hahm;
            //   }
            // }
          } // loop over ha

          // calculate number on ART at end of ts, based on number or percent
          Type artnum_hts = 0.0;
          if (dt*(hts+1) < 0.5) {
            if ( (!art15plus_isperc(g, t-2)) & (!art15plus_isperc(g, t-1)) ){ // both numbers
              artnum_hts = (0.5-dt*(hts+1))*art15plus_num(g, t-2) + (dt*(hts+1)+0.5)*art15plus_num(g, t-1);
            } else if (art15plus_isperc(g, t-2) & art15plus_isperc(g, t-1)){ // both percentages
              Type artcov_hts = (0.5-dt*(hts+1))*art15plus_num(g, t-2) + (dt*(hts+1)+0.5)*art15plus_num(g, t-1);
              artnum_hts = artcov_hts * (Xart_15plus + Xartelig_15plus);
            } else if ( (!art15plus_isperc(g, t-2)) & art15plus_isperc(g, t-1)) { // transition from number to percentage
              Type curr_coverage = Xart_15plus / (Xart_15plus + Xartelig_15plus);
              Type artcov_hts = curr_coverage + (art15plus_num(g, t-1) - curr_coverage) * dt / (0.5-dt*hts);
              artnum_hts = artcov_hts * (Xart_15plus + Xartelig_15plus);
            }
          } else {
            if( (!art15plus_isperc(g, t-1)) & (!art15plus_isperc(g, t)) ){ // both numbers
              artnum_hts = (1.5-dt*(hts+1))*art15plus_num(g, t-1) + (dt*(hts+1)-0.5)*art15plus_num(g, t);
            } else if(art15plus_isperc(g, t-1) & art15plus_isperc(g, t)){ // both percentages
              Type artcov_hts = (1.5-dt*(hts+1))*art15plus_num(g, t-1) + (dt*(hts+1)-0.5)*art15plus_num(g, t);
              artnum_hts = artcov_hts * (Xart_15plus + Xartelig_15plus);
            } else if( (!art15plus_isperc(g, t-1)) & art15plus_isperc(g, t)){ // transition from number to percentage
              Type curr_coverage = Xart_15plus / (Xart_15plus + Xartelig_15plus);
              Type artcov_hts = curr_coverage + (art15plus_num(g, t) - curr_coverage) * dt / (1.5-dt*hts);
              artnum_hts = artcov_hts * (Xart_15plus + Xartelig_15plus);
            }
          }

	  // Desired number to initiate on ART
          Type artinit_hts = artnum_hts > Xart_15plus ? artnum_hts - Xart_15plus : 0.0;

          // Use mixture of eligibility and expected mortality for initiation distribution

          for(int ha = hIDX_15PLUS; ha < hAG; ha++) {
	    for(int hm = anyelig_idx; hm < hDS; hm++) {

              if (Xartelig_15plus > 0.0) {
                Type artinit_hahm = artinit_hts * artelig_hahm(hm, ha-hIDX_15PLUS) * ((1.0 - art_alloc_mxweight)/Xartelig_15plus + art_alloc_mxweight * cd4_mort(hm, ha, g) / expect_mort_artelig15plus);
                if (artinit_hahm > artelig_hahm(hm, ha-hIDX_15PLUS)) {
                  artinit_hahm = artelig_hahm(hm, ha-hIDX_15PLUS);
                }
                if (artinit_hahm > hivstrat_adult(hm, ha, g, t) + dt * grad(hm, ha, g)) {
                  artinit_hahm = hivstrat_adult(hm, ha, g, t) + dt * grad(hm, ha, g);
                }
                grad(hm, ha, g) -= artinit_hahm / dt;
		gradART(ART0MOS, hm, ha, g) += artinit_hahm / dt;
                artinit(hm, ha, g, t) += artinit_hahm;
	      }
            }
          }
        } // end ART initiation

        for(int g = 0; g < NG; g++) {
          for(int ha = 0; ha < hAG; ha++) {
            for(int hm = everARTelig_idx; hm < hDS; hm++) {
              for(int hu = 0; hu < hTS; hu++) {
                artstrat_adult(hu, hm, ha, g, t) += dt * gradART(hu, hm, ha, g);
              }
            }
          }
        }



      } // if(t >= t_ART_start)


      for(int g = 0; g < NG; g++) {
        for(int ha = 0; ha < hAG; ha++) {
          for(int hm = 0; hm < hDS; hm++) {
            hivstrat_adult(hm, ha, g, t) += dt * grad(hm, ha, g);

          }
        }
      }

      // remove hivdeaths from hivpop and totpop
      for(int g = 0; g < NG; g++){

        // sum HIV+ population size in each hivpop age group
        TensorFixedSize<Type, Sizes<hAG>> hivpop_ha;
        int a = pIDX_HIVADULT;
        for(int ha = 0; ha < hAG; ha++){
          hivpop_ha(ha) = 0.0;
          for(int i = 0; i < hAG_SPAN[ha]; i++){
            hivpop_ha(ha) += hivpop1(a, g, t);
            a++;
          }
        }

        // remove hivdeaths proportionally to age-distribution within each age group
        a = pIDX_HIVADULT;
        for(int ha = 0; ha < hAG; ha++){
          if(hivpop_ha(ha) > 0){
            Type hivqx_ha = hivdeaths_ha(ha, g) / hivpop_ha(ha);
            for(int i = 0; i < hAG_SPAN[ha]; i++){
              Type hivdeaths_a = hivpop1(a, g, t) * hivqx_ha;
              hivdeaths(a, g, t) += hivdeaths_a;
              totpop1(a, g, t) -= hivdeaths_a;
              hivpop1(a, g, t) -= hivdeaths_a;
              a++;
            }
          } else {
            a += hAG_SPAN[ha];
          }  // end if(pop_ha[ha] > 0)
        }
      }



    }

    //ROB: population adjustment fx (start)
    // loop hiv_steps_per_year
    // adjust population to match target population size
    // TensorFixedSize<Type, Sizes<hAG, NG>> popadjprob;
    // TensorFixedSize<Type, Sizes<hAG, NG>> hivpopadjprob;
    for(int g = 0; g < NG; g++){
      for(int ha = 0; ha < (pAG-1); ha++){
        //popadjprob(ha, g) = 0.0;
        //popadjprob(ha, g) = basepop(ha, g, t) / totpop1(ha, g, t);
        // maybe need to change this for coarse age groups
        // hivpopadjprob(ha, g) = popadjprob(ha, g) ;


     hivpop1(ha, g, t) =  hivpop1(ha, g, t) * basepop(ha, g, t) / totpop1(ha, g, t);
     for(int hm = 0; hm < hDS; hm++){
       hivstrat_adult(hm, ha, g, t) = hivstrat_adult(hm, ha, g, t) * basepop(ha, g, t) / totpop1(ha, g, t);
     }
    totpop1(ha, g, t) = basepop(ha, g, t);

        //hivpop1(ha, g, t) = hivpopadjprob(ha, g) * hivpop1(ha, g, t);
        //if (t >= t_ART_start) {
        //to do: need to add in scalar to ART population
        //}

      }
    }

    //going to put ART here for right now, it might be the wrong place eventually


   for(int g = 0; g < NG; g++){
      for(int a = 0; a < pAG; a++){
        //changing this because eventually the tot pop will be basepop
       hivnpop1(a, g, t) = (totpop1(a, g, t)) - hivpop1(a, g, t);

      }
    }
   //ROB: population adjustment fx (end)



   //ROB: births to HIV+ women (start)
   double asfr_sum;
   asfr_sum = 0.0 ;
   for(int af = 0; af < 35; af++) {
     asfr_sum += asfr(af, t);
   }

   double births_sum;
   births_sum = births(t);

   double birthsCurrAge;
   double birthsHE;
   birthsHE = 0.0;
   double birthsHE_15_24;
   birthsHE_15_24 = 0.0;

   for(int af = 0; af < 35; af++) {
       double nHIVcurr;
       double nHIVlast;
       double totpop;
   nHIVcurr = 0.0;
   nHIVlast = 0.0;
     for(int hm = 0; hm < hDS; hm++){
       nHIVcurr += hivstrat_adult(hm, af, 1, t);
       nHIVlast += hivstrat_adult(hm, af, 1, t-1);
       for(int hu = 0; hu < hTS; hu++){
         nHIVcurr += artstrat_adult(hu, hm, af, 1, t);
         nHIVlast += artstrat_adult(hu, hm, af, 1, t-1);
       }
     }

   totpop = nHIVcurr + hivnpop1(af + 15,1,t);

   double prev;
   prev = nHIVcurr / totpop;

   double df;
   df = 0.0;
   double laf;
   laf = 1;
     for(int hm = 0; hm < hDS; hm++){
          df += laf * fert_mult_by_age(af) * fert_mult_offart(hm) * ((hivstrat_adult(hm, af, 1, t) + hivstrat_adult(hm, af, 1, t-1)) / 2);
       //women on ART less than 6 months use the off art fertility multiplier

          df += laf * fert_mult_by_age(af) * fert_mult_offart(hm) * ((artstrat_adult(0, hm, af, 1, t) + artstrat_adult(0, hm, af, 1, t-1)) / 2);
       for(int hu = 1; hu < hTS; hu++){
          df += laf * fert_mult_onart(af) * ((artstrat_adult(hu, hm, af, 1, t) + artstrat_adult(hu, hm, af, 1, t-1)) / 2);
       }
     }


   if(nHIVcurr > 0){
     df = df / ((nHIVcurr + nHIVlast) / 2);
   }else{
     df = 1;
   }


    birthsCurrAge = (nHIVcurr + nHIVlast) / 2 * tfr(t) * df / (df * prev + 1 - prev) *  asfr(af, t) / asfr_sum ;
     birthsHE += birthsCurrAge;
     if(af < 9){
       birthsHE_15_24 += birthsCurrAge;
     }



   }

  // birthsHE = std::round(birthsHE * 100000.0) / 100000.0;
   hiv_births(t) = birthsHE;
   //ROB: births to HIV+ women (end)

   //ROB: PMTCT need  (start)
//TODO ABORTION
   needPMTCT(t) = birthsHE;
   HIVPregwomen(t) = birthsHE_15_24;



   double sumARV;
   sumARV = 0.0;
   //just doing this for percentages rn
   for(int hp = 0; hp < hPS; hp++){
    sumARV += pmtct(hp,t,1);
   }
   double numPMTCT;
   numPMTCT = sumARV ;// + PATIENTS REALLOCATED

   //Maggie: TODO: make pmtct coverages variables that are automatically converted to % covered, I think that should work?
   //Cap so that it snaps the distribution to one if needed, with the distribution matching the input numbers
   TensorFixedSize<Type, Sizes<hPS>> pmtct_cov;
   for(int hp = 0; hp < hPS; hp++){
     if(pmtct_input_isperc(t)){
       pmtct_cov(hp) = pmtct(hp, t, 1);
     }else{
       pmtct_cov(hp) = pmtct(hp, t, 0) / sumARV;
     }
   }


   //need to add in reallocated patients as well
   if(pmtct_input_isperc(t)){ //input is percent
     OnPMTCT(t) = (needPMTCT(t) * sumARV) > 0 ? needPMTCT(t) * sumARV : 0.0;
    }else{//input is number
     OnPMTCT(t) = sumARV;
    }

   double need;
   need = needPMTCT(t) > numPMTCT ? needPMTCT(t) : numPMTCT;
   //ROB: PMTCT need  (end)


   double PMTCT_NONE;
   PMTCT_NONE = (1 - sumARV) > 0 ? 1 - sumARV : 0;
   //TOOO: ART dropout and then recalculate how many women aren't on any PMTCT

   //ROB: calcBF (start)
   //Proportion of pregnant women by CD4 count
   double sum1;
   double sum2;
   double sum3;
   sum1 = 0.0;
   sum2 = 0.0;
   sum3 = 0.0;

   for(int af = 0; af < 35; af++) {

    // for(int hm = 0; hm < hDS; hm++){
        sum3 += hivstrat_adult(4,af,1,t) + hivstrat_adult(5,af,1,t) + hivstrat_adult(6,af,1,t) ;
        sum1 += hivstrat_adult(3,af,1,t) + hivstrat_adult(2,af,1,t) ;
        sum2 += hivstrat_adult(0,af,1,t) + hivstrat_adult(1,af,1,t) ;
    // }

   }

   double total;
   total = sum1 + sum2 + sum3;

   double proplt200;
   double prop200to350;
   double propgte350;



   if(total >0){
     proplt200 = sum3/ total;
     prop200to350 = sum1 / total;
     propgte350 = sum2 / total;
   }else{
     proplt200 = 0;
     prop200to350 = 1;
     propgte350 = 0;
   }

   double proplte350 ;
   proplte350 = proplt200 + prop200to350;
   //adjust option A and option B
   //Option A and B were only authorized for women with greater than 350 CD4, so if the percentage of women
   //on option A/B > the proportion of women in this cd4 category, we assume that some must have a cd4 less than 350
   //option AB will be less effective for these women so we adjust for that
   double excessratio;
   double optA_tr;
   double optB_tr;
   //pmtct(2,t, 1) = Option A and pmtct(3, 4, 1) = option B
   if((pmtct(0,t,1) + pmtct(1,t,1)) > propgte350){
     if(propgte350 > 0){
       excessratio = ((pmtct(0,t,1) + pmtct(1,t,1)) / propgte350) - 1;
     }else{
       excessratio = 0;
     }
     //pmtct_mtct(0,0,0) option A, pmtct_mtct(0,1,0) is option b
     optA_tr = pmtct_mtct(0,0,0) * (1 + excessratio);
     optB_tr = pmtct_mtct(0,1,0) * (1 + excessratio);
   }

   else{
     excessratio = 0.0;
     optA_tr = pmtct_mtct(0,0,0) * (1 + excessratio);
     optB_tr = pmtct_mtct(0,1,0) * (1 + excessratio);
   }



   //Calculate transmission rate
   //this should be aligned or have some more strategic indexing later on
   double ptr1;
   ptr1 = 0;
   //NOTE: DROPOUT NOT YET WORKING
   double onARTretained;
   double startingARTretained;
  // onARTretained = pmtct(4,t,1) * (1 - (1 - (pmtct_dropout(0,t) /100) )* 5/12);
  // startingARTretained = pmtct(5,t,1) * (1 - (1 - (pmtct_dropout(1,t) /100) )* 5/12);
   onARTretained = pmtct(4,t,1) * ((pmtct_dropout(0,t) /100));
   startingARTretained = pmtct(5,t,1) * ((pmtct_dropout(1,t) /100));
   ptr1 = pmtct(3-1,t,1) * 0.075 +    //sdnvp
     pmtct(4-1,t,1) *  0.022 + // dual arv
     pmtct(1-1,t,1) * optA_tr + pmtct(2-1,t,1) *  optB_tr  +
     onARTretained * 0.0026 +
     startingARTretained * 0.014 +
     pmtct(6,t,1) *  0.082;

   double percent_in_program;
   percent_in_program = pmtct(0,t,1) + pmtct(1,t,1) + pmtct(2,t,1) + pmtct(3,t,1) + onARTretained + startingARTretained + pmtct(6,t,1);

   double ptr2;
   if(percent_in_program > 0){
     ptr2 = ptr1 / percent_in_program;
   }else{
     ptr2 = 0;
   }
   PMTCT_NONE = 1 - percent_in_program;

double proplt200_tr =  0.37;
double prop200to350_tr = 0.27;
double propgte350_tr = 0.15;
   //Add in women not receiving any form of prophylaxis
   if(total > 0){
     ptr1 = ptr1 + PMTCT_NONE * (proplt200 * proplt200_tr + prop200to350 * prop200to350_tr + propgte350 * propgte350_tr);
   }else{
     ptr1 = ptr1 ;
   }
   double ptr3;
   ptr3 = ptr1;

   //Add in transmission due to incident infections
   sum2 = 0.0; //HIV negative 15-49 women weighted for ASFR
   sum3 = 0.0; //newly infected 15-49 women, weighted for ASFR

   for(int af = 0; af < 35; af++) {
     sum2 += asfr(af, t) / asfr_sum  * hivnpop1(af + 15,1,t) ;
     sum3 +=  asfr(af, t) / asfr_sum  * infections(af + 15,1,t) ;
   }

   double IncRate ;
  // if(sum2 > need){
     IncRate = sum3 / sum2;
 //  }else{
  //   IncRate = 0;
  // }

   double v3;
   //need to pull incidence infection mtct into the input object
   v3 = IncRate * (9/12) * (births_sum - need) *0;  //.181;

   if(need > 0){
     ptr1 = ptr1 + v3 / need;
   }
   //CHANGED SOMETHING HERE
   double hivpos_births;

   hivpos_births = birthsHE * ptr1;


   double existinghivbirths;
   existinghivbirths = birthsHE * ptr3;

   double birthsHE_bf;
   birthsHE_bf = birthsHE - hivpos_births;



    //BREASTFEEDING TRANSMISSION
   double v2;
   v2 = 0.0;
   //insert hBF as breastfeeding durations
   //bf_duration is % of women no longer breast feeding
   for(int bf = 0; bf < hBF; bf++){
     v2 += IncRate / 12 * 2 * (1 - bf_duration(bf, t, 0));
   }
   double v4;
  // v4 = v2 * 0.269;
   v4 = 0.0;
   //Incident infections are hiv+mothers minus hiv births * v4 (which has already been adjusted for prevalence)
   double IncidentInfectionsBF;
   IncidentInfectionsBF = (births_sum - needPMTCT(t)) * v4;


   //baseline bftr = 0
   double bftr;
   bftr = 0.0;

   double optA_trbf;
   double optB_trbf;
   double excess;
   optA_trbf = 0;
   optB_trbf = 0;
   excess = 0;

   if(propgte350 > 0){
     if((pmtct(0,t,1) + pmtct(1,t,1)) > propgte350){
    //if((pmtct(0,t,1) + pmtct(1,t,1) - pmtct(4,t,1) - pmtct(5,t,1) - pmtct(6,t,1)) > propgte350){
      // excess = pmtct(0,t,1) + pmtct(1,t,1) - pmtct(4,t,1) - pmtct(5,t,1) - pmtct(6,t,1) - propgte350;
       excess = pmtct(0,t,1) + pmtct(1,t,1) - propgte350;
       optA_trbf = (propgte350 * pmtct_mtct(4,0,1)) + excess * (1.45 / 0.46) * pmtct_mtct(4,0,1) / (propgte350 + excess);
       optB_trbf = (propgte350 * pmtct_mtct(4,1,1)) + excess * (1.45 / 0.46) * pmtct_mtct(4,1,1) / (propgte350 + excess);

     }else{
       optA_trbf = pmtct_mtct(4,0,1);
       optB_trbf = pmtct_mtct(4,1,1);

     }
   }else{
     optA_trbf = pmtct_mtct(4,0,1);
     optB_trbf = pmtct_mtct(4,1,1);
   }
   //bftr from birth to <6 months
   double trt_pct;


   double bftr_1;
   double NoPMTCT_bf;
   double NewInfBFLt6;
   NewInfBFLt6 = 0.0;
   bftr_1 = 0.0;


   double tr_lt350_bf = 0.0089;
   double tr_gt350_bf = 0.0081;

   //Note that incidence isn't dependent
   for(int bf = 0; bf < 3; bf++){
     //ptr3 is the transmission that has already occurred due to perinatal transmission
     //NoPMTCT_bf is the percentage of women who are still vulnerable to HIV transmission to their babies
     NoPMTCT_bf = 1 - ptr3 - bftr_1;

      for(int hp = 0; hp < 7; hp++){
        //hp = 0 is option A
        //Dropout not used for option A
        if(hp == 0){
          trt_pct = optA_trbf * pmtct(0,t,1);// * (1 - (pmtct_dropout(2,t) / 100) * 2);
          bftr_1 += trt_pct * 2 * (1 - bf_duration(bf, t, 1)) ;
          NoPMTCT_bf -= pmtct(0,t,1);
        }
        //hp = 1 is option B
        //Dropout not used for option B
        if(hp == 1){
          trt_pct = optB_trbf * pmtct(1,t,1) ; //* (1 - (pmtct_dropout(3,t) / 100) * 2);
          bftr_1 +=  trt_pct * 2 * (1 - bf_duration(bf, t, 1)) ;
          NoPMTCT_bf -=  pmtct(1,t,1) ;
        }
        if(hp > 3){
          //This isn't the correct implementation of dropout, see page 174 of Spectrum manual
         // trt_pct = pmtct(hp,t,1) * 1 / exp((bf+1) * 2 * log(1 + pmtct_dropout(4,t) / 100)) ;
          if(bf > 0){
            trt_pct = pmtct(hp,t,1) * (pow(1 - pmtct_dropout(4,t) / 100 * 2, bf))  ;
          }else{
            trt_pct = pmtct(hp,t,1);
          }
          bftr_1 += trt_pct * 2 * (1 - bf_duration(bf, t, 1)) * pmtct_mtct(4,hp,1);
          NoPMTCT_bf -= trt_pct ;

        }
      }
         if(NoPMTCT_bf < 0){
           NoPMTCT_bf = 0;
         }


         //No treatment
         if(bf_duration(bf, t, 0) < 1){
           if(total > 0){
             bftr_1 +=  NoPMTCT_bf * (1 - bf_duration(bf, t, 0)) * (2 * (1 - propgte350) * tr_lt350_bf + 2 * propgte350 * tr_gt350_bf);
           }
         }


      if(bf < 1){
         bftr_1 = bftr_1/ 4;
       }


   }
   NewInfBFLt6 = birthsHE  * bftr_1;



   //bftr from 6-12 months
   double bftr_2;
   double NewInfBFgte6;
   NewInfBFgte6 = 0.0;
   bftr_2 = 0.0;
   for(int bf = 3; bf < 6; bf++){

       NoPMTCT_bf = 1 - ptr3 - bftr_2;


     for(int hp = 0; hp < 7; hp++){
       //hp = 0 is option A
       //Dropout not used for option A
       if(hp == 0){
         trt_pct = optA_trbf * pmtct(0,t,1) ;//* (1 - (pmtct_dropout(2,t) / 100) * 2);
         bftr_2 +=  trt_pct * 2 * (1 - bf_duration(bf, t, 1))  ;
         NoPMTCT_bf -=  pmtct(0,t,1);
       }
       //hp = 1 is option B
       //Dropout not used for option B
       if(hp == 1){
         trt_pct =  optB_trbf * pmtct(1,t,1)   ;//* (1 - (pmtct_dropout(3,t) / 100) * 2);
         bftr_2 +=  trt_pct * 2 * (1 - bf_duration(bf, t, 1)) ;
         NoPMTCT_bf -= pmtct(1,t,1);
       }
       if(hp > 3){
         trt_pct = pmtct(hp,t,1) * (pow(1 - pmtct_dropout(4,t) / 100 * 2, (bf))) ;
         bftr_2 +=  trt_pct * 2 * (1 - bf_duration(bf, t, 1)) * pmtct_mtct(4,hp,1);
         NoPMTCT_bf -=  trt_pct;
       }
     }

       if(NoPMTCT_bf < 0){
         NoPMTCT_bf = 0;
       }
       //No treatment
       bftr_2 +=  NoPMTCT_bf * (1 - bf_duration(bf, t, 0)) *  (2 * (1 - propgte350) * tr_lt350_bf + 2 * propgte350 * tr_gt350_bf);

   }
   NewInfBFgte6 = birthsHE * bftr_2;
//   NewInfBFgte6 = std::round(NewInfBFgte6 * 100000.0) / 100000.0;


   //bftr from 12+months
   double NewInfBFgte12;
   double bftr_3;
   NewInfBFgte12 = 0.0;

   bftr_3 = 0.0;

   for(int bf = 6; bf < 12; bf++){

       NoPMTCT_bf = 1 - ptr3 - bftr_3;

     for(int hp = 0; hp < 7; hp++){
       //hp = 0 is option A
       //Dropout not used for option A
       if(hp == 0){
         trt_pct = optA_trbf * pmtct(0,t,1)  ;//* (1 - (pmtct_dropout(2,t) / 100) * 2);
         bftr_3 +=  trt_pct * 2 * (1 - bf_duration(bf, t, 1))  ;
         NoPMTCT_bf -=  pmtct(0,t,1);
       }
       //hp = 1 is option B
       //Dropout not used for option B
       if(hp == 1){
         trt_pct =  optB_trbf * pmtct(1,t,1) ;// * (1 - (pmtct_dropout(3,t) / 100) * 2 );
         bftr_3 +=  trt_pct * 2 * (1 - bf_duration(bf, t, 1));
         NoPMTCT_bf -= pmtct(1,t,1);
       }
       if(hp > 3){
         trt_pct = pmtct(hp,t,1)  *  pow(1 - pmtct_dropout(4,t) / 100 * 2, 5) *  pow(1 - pmtct_dropout(5,t) / 100 * 2, (bf-5)) ;
         bftr_3 +=  trt_pct * 2 * (1 - bf_duration(bf, t, 1)) * pmtct_mtct(4,hp,1);
         NoPMTCT_bf -=  trt_pct;
       }
     }
     if(NoPMTCT_bf < 0){
       NoPMTCT_bf = 0;
     }
     //No treatment
       bftr_3 +=   NoPMTCT_bf * (1 - bf_duration(bf, t, 0)) *  (2 * (1 - propgte350) * tr_lt350_bf + 2 * propgte350 * tr_gt350_bf);
   }


   //bftr from 24-36
   double NewInfBFgte24;
   double bftr_4;
   bftr_4 = 0.0;
   NewInfBFgte24 = 0.0;
   for(int bf = 12; bf < hBF; bf++){
    NoPMTCT_bf = 1 - ptr3 - bftr_4;

     for(int hp = 0; hp < 7; hp++){
       //hp = 0 is option A
       //Dropout not used for option A
       if(hp == 0){
         trt_pct = optA_trbf * pmtct(0,t,1) ;//* (1 - (pmtct_dropout(2,t) / 100) * 2);
         bftr_4 +=  trt_pct * 2 * (1 - bf_duration(bf, t, 1)) ;
         NoPMTCT_bf -=  pmtct(0,t,1);
       }
       //hp = 1 is option B
       //Dropout not used for option B
       if(hp == 1){
         trt_pct = optB_trbf * pmtct(1,t,1)  ;// * (1 - (pmtct_dropout(3,t) / 100) * 2);
         bftr_4 +=  trt_pct * 2 * (1 - bf_duration(bf, t, 1)) ;
         NoPMTCT_bf -= pmtct(1,t,1);
       }
       if(hp > 3){
         trt_pct = pmtct(hp,t,1)  * pow(1 - pmtct_dropout(4,t) / 100 * 2, 5) *  (pow(1 - pmtct_dropout(5,t) / 100 * 2, (bf-5)));
         bftr_4 +=  trt_pct * 2 * (1 - bf_duration(bf, t, 1)) * pmtct_mtct(4,hp,1);
         NoPMTCT_bf -=  trt_pct;
       }
     }
     if(NoPMTCT_bf < 0){
       NoPMTCT_bf = 0;
     }
     //No treatment
     bftr_4 +=  NoPMTCT_bf * (1 - bf_duration(bf, t, 0)) *  (2 * (1 - propgte350) * tr_lt350_bf + 2 * propgte350 * tr_gt350_bf);
   }
   //ROB: calcBF (end)




   //ROB: distribute MTCT infections (start)
   //will need to add in nosocomial
   TensorFixedSize<Type, Sizes<5, NG>> temp_inf;
   temp_inf.setZero();

   for(int g = 0; g < NG; g++){
     //not sure if hiv free surv is needed here
     temp_inf(0, g) += hivpos_births * births_sex_prop(g, t);
     //temp_inf(1, g) = (NewInfBFLt6) * births_sex_prop(g, t);
     temp_inf(1, g) = (NewInfBFLt6 + IncidentInfectionsBF) * births_sex_prop(g, t);
     temp_inf(2, g) = (NewInfBFgte6) * births_sex_prop(g, t);

     hivnpop1(0,g,t) -= temp_inf(0, g);
     hivnpop1(0,g,t) -= temp_inf(1, g);


     double nNeg;
     nNeg = hivnpop1(1,0,t) + hivnpop1(1,1,t) ;
     NewInfBFgte12 = 0.0;
     if(nNeg > 0 ){
       NewInfBFgte12 = birthsHE * bftr_3 * (hivnpop1(1,g,t) / nNeg);
     }
     //Remove those who were infected
     hivnpop1(1,g,t) -= NewInfBFgte12;
     temp_inf(3, g) = NewInfBFgte12 ;

     infections(0,g,t) += temp_inf(0,g) + temp_inf(1,g) + temp_inf(2,g);
     infections(1,g,t) += temp_inf(3,g);
     NewInfBFgte24 = 0.0;
     nNeg = hivnpop1(2,0,t) + hivnpop1(2,1,t) ;
     if(nNeg > 0 ){
       NewInfBFgte24 = birthsHE * bftr_4 * (hivnpop1(2,g,t) / nNeg);
     }
     //Remove those who were infected
     hivnpop1(2,g,t) -= NewInfBFgte24;

     temp_inf(4, g) = NewInfBFgte24;
     infections(2,g,t) += temp_inf(4, g);

   }

   //ROB: nosocomial infections (start)
   //distribute across eligible ages, right now just going to hardcode
   for(int g = 0; g < NG; g++){
     for(int af = 0; af < 5; af++){
       if(paed_incid_input(t) > 0){
         infections(af, g, t) = 0 ;
         infections(af, g, t) = paed_incid_input(t) / 10;
         hivpop1(af, g, t) += infections(af, g, t);

         for(int hm = 0; hm < hDS; hm++){
           for(int cat = 0; cat < 1; cat++){
             //putting them all in perinatal hTM to match spec nosocomial
             hivstrat_paeds(hm, cat, af, g, t) += paed_cd4_dist(hm) > 0 ? infections(af, g, t) * paed_cd4_dist(hm) : 0.0;
           }
         }
       }
     }

     for(int hm = 0; hm < hDS; hm++){
       for(int cat = 0; cat < 1; cat++){
         //   hivstrat_paeds(hm, cat, 0, g, t) += paed_cd4_dist(hm) > 0 ? infections(0, g, t) * paed_cd4_dist(hm) : 0.0;

       }
     }
   }
   //ROB: nosocomial infections (end)

   for(int hm = 0; hm < hDS; hm++){
     for(int g = 0; g < NG; g++){
       hivstrat_paeds(hm, 0, 0, g, t) +=  temp_inf(0,g)  * paed_cd4_dist(hm) ;
       hivstrat_paeds(hm, 1, 0, g, t) +=  temp_inf(1,g) * paed_cd4_dist(hm) ;
       hivstrat_paeds(hm, 2, 0, g, t) +=  temp_inf(2,g) * paed_cd4_dist(hm) ;
       hivstrat_paeds(hm, 3, 1, g, t) +=  temp_inf(3,g) * paed_cd4_dist(hm) ;
       hivstrat_paeds(hm, 3, 2, g, t) +=  temp_inf(4, g) * paed_cd4_dist(hm) ;
     }
   }
   //ROB: distribute MTCT infections (end)


   //ROB: paediatric natural history (start)
   for(int g = 0; g < NG; g++){
      for(int hm = 0; hm < hDS; hm++){
        for(int af = 0; af < 5; af++){
          for(int cat = 0; cat < hTM; cat++){
            deaths_paeds(hm, cat, af, g, t) += hivstrat_paeds(hm, cat, af, g, t) - (1 - ctx_effect * ctx_val(t)) * hivstrat_paeds(hm, cat, af, g, t) * paed_cd4_mort(hm, cat, af) ;
            aidsdeaths_noart_paed(hm, cat, af, g, t) += (1 - ctx_effect * ctx_val(t)) * hivstrat_paeds(hm, cat, af, g, t) * paed_cd4_mort(hm, cat, af)  ;
          }
        }
      }
    }



    for(int g = 0; g < NG; g++){
      for(int hm = 0; hm < hDS_adol; hm++){
        for(int af = 5; af < pIDX_HIVADULT; af++){
          for(int cat = 0; cat < hTM; cat++){
             deaths_paeds(hm, cat, af, g, t) += hivstrat_paeds(hm, cat, af, g, t) - (1 - ctx_effect * ctx_val(t)) * hivstrat_paeds(hm, cat, af, g, t) * adol_cd4_mort(hm, cat, af - 5);
             aidsdeaths_noart_paed(hm, cat, af, g, t) +=  (1 - ctx_effect * ctx_val(t)) * hivstrat_paeds(hm, cat, af, g, t) * adol_cd4_mort(hm, cat, af - 5); // output hiv deaths, aggregated across transmission category
          }
        }
      }
    }

    //progress through CD4 categories
    for(int g = 0; g < NG; g++){
      for(int hm = 1; hm < hDS; hm++){
        for(int af = 0; af < 5; af++){
          for(int cat = 0; cat < hTM; cat++){
            // grad_paeds(hm - 1, cat, af, g, t) -= ((hivstrat_paeds(hm-1, cat, af-1, g, t-1) + artstrat_paeds(hm-1, cat, af-1, g, t-1)) > 0) ? (deaths_paeds(hm - 1, cat, af, g, t) * paed_cd4_prog(hm - 1) + hivstrat_paeds(hm - 1, cat, af, g, t) * paed_cd4_prog(hm - 1)) / 2: 0.0; //moving to next cd4 category
            // grad_paeds(hm, cat, af, g, t) += ((hivstrat_paeds(hm-1, cat, af-1, g, t-1) + artstrat_paeds(hm-1, cat, af-1, g, t-1)) > 0) ? (deaths_paeds(hm - 1, cat, af, g, t) * paed_cd4_prog(hm - 1) + hivstrat_paeds(hm - 1, cat, af, g, t) * paed_cd4_prog(hm - 1)) / 2 : 0.0; //moving into this cd4 category
           // if(paed_cd4_mort(hm, cat, af) > 0 or paed_cd4_mort(hm-1, cat, af) > 0){
              grad_paeds(hm - 1, cat, af, g, t) -=  (deaths_paeds(hm - 1, cat, af, g, t) * paed_cd4_prog(hm - 1) + hivstrat_paeds(hm - 1, cat, af, g, t) * paed_cd4_prog(hm - 1)) / 2; //moving to next cd4 category
              grad_paeds(hm, cat, af, g, t) += (deaths_paeds(hm - 1, cat, af, g, t) * paed_cd4_prog(hm - 1) + hivstrat_paeds(hm - 1, cat, af, g, t) * paed_cd4_prog(hm - 1)) / 2; //moving into this cd4 category
           // }
          }
        }
      }
    }

    //progress through CD4 categories
    for(int g = 0; g < NG; g++){
      for(int hm = 1; hm < hDS_adol; hm++){
        for(int af = 5; af < pIDX_FERT; af++){
          for(int cat = 0; cat < hTM; cat++){
        //   if(adol_cd4_mort(hm, cat, af - 5) > 0 or adol_cd4_mort(hm - 1, cat, af - 5) > 0){
              grad_paeds(hm - 1, cat, af, g, t) -= (deaths_paeds(hm - 1, cat, af, g, t) * adol_cd4_prog(hm - 1) + hivstrat_paeds(hm - 1, cat, af, g, t) * adol_cd4_prog(hm - 1)) / 2; //moving to next cd4 category
              grad_paeds(hm, cat, af, g, t) += (deaths_paeds(hm - 1, cat, af, g, t) * adol_cd4_prog(hm - 1) + hivstrat_paeds(hm - 1, cat, af, g, t) * adol_cd4_prog(hm - 1)) / 2; //moving into this cd4 category
        //    }
          }
        }
      }
    }
    //ROB: paediatric natural history (end)


    //ROB: paediatric offART HIV mortality (start)
    //HIV RELATED MORTALTITY
    for(int g = 0; g < NG; g++){
      for(int hm = 0; hm < hDS; hm++){
        for(int af = 0; af < 5; af++){
          for(int cat = 0; cat < hTM; cat++){
            grad_paeds(hm, cat, af, g, t) -= (1 - ctx_effect * ctx_val(t)) * hivstrat_paeds(hm, cat, af, g, t) * paed_cd4_mort(hm, cat, af)  ;
          }
        }
      }
    }

    //not accounting for ctx 5 year
    for(int g = 0; g < NG; g++){
      for(int hm = 0; hm < hDS_adol; hm++){
        for(int af = 5; af < pIDX_FERT; af++){
          for(int cat = 0; cat < hTM; cat++){
            //TO DO: add in cotrim effects
            grad_paeds(hm, cat, af, g, t) -= (1 - ctx_effect * ctx_val(t)) * hivstrat_paeds(hm, cat, af, g, t) * adol_cd4_mort(hm, cat, af - 5);
          }
        }
      }
    }
    //ROB: paediatric offART HIV mortality (end)

  //not here mkw

    for(int g = 0; g < NG; g++){
      for(int hm = 0; hm < hDS; hm++){
        for(int af = 0; af < pIDX_FERT; af++){
          for(int cat = 0; cat < hTM; cat++){
            hivstrat_paeds(hm, cat, af, g, t) += grad_paeds(hm, cat, af, g, t) ;
          }
        }
      }
    }

//here mkw

  //ROB: paediatric ART retention (start)
    //LTFU, note this may be in the wrong place. Also assuming that distribution will be according to the distribution
    //of transmission types in the current HIV population, that might be wrong
    //Note, I've verified that this iswrong
   TensorFixedSize<Type, Sizes<hTM, pIDX_HIVADULT, NG>> tothiv_pop_bycat;
   for(int hm = 0; hm < hDS; hm++){
     for(int af = 0; af < pIDX_FERT; af++){
       for(int g = 0; g < NG; g++){
         for(int cat = 0; cat < hTM; cat++){
           tothiv_pop_bycat(cat, af, g) += hivstrat_paeds(hm, cat, af, g, t);
         }
       }
     }
   }

  TensorFixedSize<Type, Sizes<pIDX_HIVADULT, NG>> tothiv_pop;
     for(int af = 0; af < pIDX_FERT; af++){
       for(int g = 0; g < NG; g++){
         for(int cat = 0; cat < hTM; cat++){
           tothiv_pop(af, g) += tothiv_pop_bycat(cat, af, g);

         }
       }
     }

  TensorFixedSize<Type, Sizes<hDS, pIDX_HIVADULT, NG>> art_total;
   for(int hm = 0; hm < hDS; hm++){
     for(int af = 0; af < pIDX_FERT; af++){
       for(int g = 0; g < NG; g++){
         for(int cat = 0; cat < hTM; cat++){
           for(int hu = 0 ; hu < hTS; hu ++){
             art_total(hm, af, g) += artstrat_paeds(hu, hm, af, g, t);
           }
           if(paed_art_ltfu(t) > 0 & tothiv_pop(af,g) > 0 ){
             hivstrat_paeds(hm, cat, af, g, t) += (tothiv_pop_bycat(cat, af, g) / tothiv_pop(af, g)) * paed_art_ltfu(t) *  art_total(hm, af, g);
           }
         }
       }
     }
   }



    for(int hm = 0; hm < hDS; hm++){
      for(int af = 0; af < pIDX_FERT; af++){
        for(int g = 0; g < NG; g++){
          for(int hu; hu < hTS; hu++){
            artstrat_paeds(hu, hm, af, g, t) -= artstrat_paeds(hu, hm, af, g, t) * paed_art_ltfu(t) ;
          }
        }
      }
    }
    //ROB: paediatric ART retention (end)


    //ROB: paediatric ART initiation (start)
    TensorFixedSize<Type, Sizes<hDS, hTM, pIDX_HIVADULT, NG>> need_art_paed;
    need_art_paed.setZero();
    //all children under a certain age eligible for ART
    for(int g = 0; g < NG; g++){
      for(int cat = 0; cat < hTM; cat++){
        for(int af = 0; af < pIDX_HIVADULT; af++){
          for(int hm = 0; hm < hDS; hm++){
            need_art_paed(hm, cat, af, g) += af < paed_art_elig_age(t) ? hivstrat_paeds(hm, cat, af, g, t) :0.0 ;
          }
        }
      }
    }


    //all children under a certain CD4 are eligible for ART, regardless of age
    for(int g = 0; g < NG; g++){
      for(int cat = 0; cat < hTM; cat++){
        for(int af = paed_art_elig_age(t); af < pIDX_HIVADULT; af++){
          for(int hm = 0; hm < hDS; hm++){
            need_art_paed(hm, cat, af, g) += (hm > (paed_art_elig_cd4(af, t) - 2)) ? hivstrat_paeds(hm, cat, af, g, t) : 0.0;
          }
        }
      }
    }


    artnum_paed(t) = 0.0;
    for(int g = 0; g < NG; g++){
        for(int af = 0; af < pIDX_HIVADULT; af++){
          for(int hm = 0; hm < hDS; hm++){
            for(int cat = 0; cat < hTM; cat++){
              artnum_paed(t) += need_art_paed(hm, cat, af, g);
            }
        }
      }
    }

    //how many should initialize ART
    for(int g = 0; g < NG; g++){
        for(int af = 0; af < pIDX_HIVADULT; af++){
          for(int hm = 0; hm < hDS; hm++){
            for(int cat = 0; cat < hTM; cat++){
              if(paed_art_val(t) > 0){
                init_art_paed(hm, cat, af, g, t) += need_art_paed(hm, cat, af, g);
                for(int dur = 0; dur < hTS; dur++){
                  init_art_paed(hm, cat, af, g, t) = init_art_paed(hm, cat, af, g, t) < 0 ? 0.0 : init_art_paed(hm, cat, af, g, t);
                }
              }
            }
        }
      }
    }

    double init_art_paed_total;
    init_art_paed_total = 0.0;
    for(int g = 0; g < NG; g++){
        for(int af = 0; af < pIDX_HIVADULT; af++){
          for(int hm = 0; hm < hDS; hm++){
            for(int cat = 0; cat < hTM; cat++){
            init_art_paed_total += init_art_paed(hm, cat, af, g, t);
          }
        }
      }
    }

  //here mkw

    for(int g = 0; g < NG; g++){
      for(int hm = 0; hm < hDS; hm++){
        for(int af = 0; af < pIDX_HIVADULT; af++){
            double death_rate ;
            death_rate = af < 5 ? mort_art_rr(0, af, t) * 0.5 * (paed_art_mort(hm, 0, af) + paed_art_mort(hm, 1, af)) : mort_art_rr(0, af, t) * 0.5 * (adol_art_mort(hm, 0, af-5) + adol_art_mort(hm, 1, af-5));
            aidsdeaths_art_paed(0,hm, af, g, t) =  death_rate * artstrat_paeds(0, hm, af, g, t)  ;
            grad_paeds_art(0,hm, af, g, t) -= aidsdeaths_art_paed(0,hm, af, g, t) ;
            artstrat_paeds(0, hm,  af, g, t) += grad_paeds_art(0, hm, af, g, t) ;
            grad_paeds_art(0, hm, af, g, t) = 0.0;
        }
      }
    }


    for(int g = 0; g < NG; g++){
      for(int hm = 0; hm < hDS; hm++){
        for(int af = 0; af < pIDX_HIVADULT; af++){
          double death_rate;
          death_rate =  af < 5 ? mort_art_rr(2, af, t) * paed_art_mort(hm, 2, af) : mort_art_rr(2, af, t) * adol_art_mort(hm, 2, af-5);
          aidsdeaths_art_paed(2,hm, af, g, t) =  artstrat_paeds(2, hm, af, g, t) * death_rate;
          grad_paeds_art(2,hm, af, g, t) -= aidsdeaths_art_paed(2,hm, af, g, t) ;
          artstrat_paeds(2, hm,  af, g, t) += grad_paeds_art(2, hm, af, g, t) ;

        }
      }
    }

    //Progress ART to the correct time on ART
    for(int hm = 0; hm < hDS; hm++){
      for(int af = 0; af < pIDX_FERT; af++){
        for(int g = 0; g < NG; g++){
          artstrat_paeds(1, hm, af, g, t) += artstrat_paeds(0, hm, af, g, t) > 0 ? artstrat_paeds(0, hm, af, g, t) : 0;
          artstrat_paeds(0, hm, af, g, t) -= artstrat_paeds(0, hm, af, g, t) > 0 ? artstrat_paeds(0, hm, af, g, t) : 0;
        }
      }
    }
   //here mkw

    if ( (!artpaeds_isperc(t)) & (!artpaeds_isperc(t-1)) ){ // both numbers
      //Remove how many that are already on ART
      artnum_paed(t) =  (paed_art_val(t) + paed_art_val(t-1)) / 2 ;
      for(int g = 0; g < NG; g++){
        for(int af = 0; af < pIDX_HIVADULT; af++){
          for(int hm = 0; hm < hDS; hm++){
            for(int dur = 0; dur < hTS; dur++){
              artnum_paed(t) -= (artstrat_paeds(dur, hm, af, g, t) )  ;

            }
          }
        }
      }

      artnum_paed(t) = init_art_paed_total < artnum_paed(t) ? init_art_paed_total : artnum_paed(t) ;
      artnum_paed(t) = artnum_paed(t) < 0 ? 0 : artnum_paed(t);

    } else if (artpaeds_isperc(t) & artpaeds_isperc(t-1)){ // both percentages
      for(int g = 0; g < NG; g++){
        for(int af = 0; af < pIDX_HIVADULT; af++){
          for(int hm = 0; hm < hDS; hm++){
            for(int dur = 0; dur < hTS; dur++){
              artnum_paed(t) += (artstrat_paeds(dur, hm, af, g, t) )  ;
              artnum_paed(t) += aidsdeaths_art_paed(dur,hm, af, g, t) ;
            }
          }
        }
      }

      artnum_paed(t) =  artnum_paed(t) * (paed_art_val(t) + paed_art_val(t-1)) / 2 ;

      //Remove how many that are already on ART
      for(int g = 0; g < NG; g++){
        for(int af = 0; af < pIDX_HIVADULT; af++){
          for(int hm = 0; hm < hDS; hm++){
            for(int dur = 0; dur < hTS; dur++){
              artnum_paed(t) -= (artstrat_paeds(dur, hm, af, g, t)) ;

            }
          }
        }
      }
      artnum_paed(t) = init_art_paed_total < artnum_paed(t) ? init_art_paed_total : artnum_paed(t) ;



    } else if (artpaeds_isperc(t) & !artpaeds_isperc(t-1)){ // num to percentage
      //Remove how many that are already on ART
      double temp ;
      temp = 0.0;
      for(int g = 0; g < NG; g++){
        for(int af = 0; af < pIDX_HIVADULT; af++){
          for(int hm = 0; hm < hDS; hm++){
            for(int dur = 0; dur < hTS; dur++){
              artnum_paed(t) += artstrat_paeds(dur, hm, af, g, t) +  aidsdeaths_art_paed(dur,hm, af, g, t)  ;

            }
          }
        }
      }
      artnum_paed(t) = (paed_art_val(t-1) + (artnum_paed(t) * paed_art_val(t))) / 2 ;

      for(int g = 0; g < NG; g++){
        for(int af = 0; af < pIDX_HIVADULT; af++){
          for(int hm = 0; hm < hDS; hm++){
            for(int dur = 0; dur < hTS; dur++){
              artnum_paed(t) -= (artstrat_paeds(dur, hm, af, g, t) )  ;


            }
          }
        }
      }

      artnum_paed(t) = artnum_paed(t) < 0 ? 0 : artnum_paed(t);
      artnum_paed(t) = init_art_paed_total < artnum_paed(t) ? init_art_paed_total : artnum_paed(t) ;

    } else if (artpaeds_isperc(t-1) & !artpaeds_isperc(t)){ //percentage to num


      for(int g = 0; g < NG; g++){
        for(int af = 0; af < pIDX_HIVADULT; af++){
          for(int hm = 0; hm < hDS; hm++){
            for(int dur = 0; dur < hTS; dur++){
              artnum_paed(t) += artstrat_paeds(dur, hm, af, g, t) ;

            }
          }
        }
      }


      artnum_paed(t) = (artnum_paed(t-1) + paed_art_val(t)) / 2 ;




      //Remove how many that are already on ART
      for(int g = 0; g < NG; g++){
        for(int af = 0; af < pIDX_HIVADULT; af++){
          for(int hm = 0; hm < hDS; hm++){
            for(int dur = 0; dur < hTS; dur++){
              artnum_paed(t) -= (artstrat_paeds(dur, hm, af, g, t))  ;

            }
          }
        }
      }


      artnum_paed(t) = artnum_paed(t) < 0 ? 0 : artnum_paed(t);
      artnum_paed(t) = init_art_paed_total < paed_art_val(t) ? init_art_paed_total : paed_art_val(t) ;


    }

    //here mkw

    double initByAge = 0.0;
    for(int g = 0; g < NG; g++){
      for(int af = 0; af < pIDX_HIVADULT; af++){
        for(int hm = 0; hm < hDS; hm++){
          for(int cat = 0; cat < hTM; cat++){
            initByAge +=  init_art_paed(hm, cat, af, g, t) * init_art_dist(af, t);
          }
        }
      }
    }
    double adj ;
    adj = initByAge == 0 ? 1 : artnum_paed(t) / initByAge ;
    for(int g = 0; g < NG; g++){
      for(int cat = 0; cat < hTM; cat++){
        for(int af = 0; af < pIDX_HIVADULT; af++){
          for(int hm = 0; hm < hDS; hm++){
            double scalar ;
            scalar = (adj * init_art_dist(af, t)) > 1 ? 1 : adj * init_art_dist(af, t);
            scalar = artnum_paed(t) > 0 ? scalar : 0.0;
            //not really sure why this needs to be here, but without it it scales down certain age groups when there is enough ART
            //scalar = (paed_art_mort(hm, cat, af) == 0) ?  1 : scalar;
            artstrat_paeds(0, hm, af, g, t) +=  scalar * init_art_paed(hm, cat, af, g, t) ;
            hivstrat_paeds(hm, cat, af, g, t) -=  scalar* init_art_paed(hm, cat, af, g, t) ;

          }
        }
      }
    }
    //ROB: paediatric ART initiation (end)

    //ROB: onART mortlity (start)
    for(int g = 0; g < NG; g++){
      for(int hm = 0; hm < hDS; hm++){
        for(int af = 0; af < pIDX_HIVADULT; af++){
          double death_rate ;
            death_rate = af < 5 ? mort_art_rr(0, af, t) * 0.5 * (paed_art_mort(hm, 0, af) + paed_art_mort(hm, 1, af)) : mort_art_rr(0, af, t) * 0.5 * (adol_art_mort(hm, 0, af-5) + adol_art_mort(hm, 1, af-5));
            aidsdeaths_art_paed(0,hm, af, g, t) =  death_rate * artstrat_paeds(0, hm, af, g, t)  ;
            grad_paeds_art(0,hm, af, g, t) -= aidsdeaths_art_paed(0,hm, af, g, t) ;
            artstrat_paeds(0, hm,  af, g, t) += grad_paeds_art(0, hm, af, g, t) ;

        }
      }
    }

    //Progress ART to the correct time on ART
    for(int hm = 0; hm < hDS; hm++){
      for(int af = 0; af < pIDX_FERT; af++){
        for(int g = 0; g < NG; g++){
            artstrat_paeds(2, hm, af, g, t) += artstrat_paeds(1, hm, af, g, t) > 0 ? artstrat_paeds(1, hm, af, g, t) : 0;
            artstrat_paeds(1, hm, af, g, t) -= artstrat_paeds(1, hm, af, g, t) > 0 ? artstrat_paeds(1, hm, af, g, t) : 0;
        }
      }
    }
    //ROB: onART mortality (end)


  // //ROB: nosocomial infections (start), this matches spectrum
  //   //distribute across eligible ages, right now just going to hardcode
  //   for(int g = 0; g < NG; g++){
  //     for(int af = 0; af < 5; af++){
  //       if(paed_incid_input(t) > 0){
  //         infections(af, g, t) = 0 ;
  //         infections(af, g, t) = paed_incid_input(t) / 10;
  //         hivpop1(af, g, t) += infections(af, g, t);
  //
  //         for(int hm = 0; hm < hDS; hm++){
  //           for(int cat = 0; cat < 1; cat++){
  //             //putting them all in perinatal hTM to match spec nosocomial
  //            hivstrat_paeds(hm, cat, af, g, t) += paed_cd4_dist(hm) > 0 ? infections(af, g, t) * paed_cd4_dist(hm) : 0.0;
  //           }
  //         }
  //       }
  //    }
  //
  //     for(int hm = 0; hm < hDS; hm++){
  //       for(int cat = 0; cat < 1; cat++){
  //       //   hivstrat_paeds(hm, cat, 0, g, t) += paed_cd4_dist(hm) > 0 ? infections(0, g, t) * paed_cd4_dist(hm) : 0.0;
  //
  //       }
  //     }
  //   }
  //   //ROB: nosocomial infections (end)

    for(int g = 0; g < NG; g++){
      for(int af =0; af < pIDX_FERT; af++){
        hivpop1(af,g,t) = 0.0;
        for(int hm = 0; hm < hDS; hm++){
          for(int cat = 0; cat < 4; cat++){
            hivpop1(af,g,t) +=  hivstrat_paeds(hm, cat, af, g, t);
          }
          for(int cat = 0; cat < hTS; cat++){
            hivpop1(af,g,t) += artstrat_paeds(cat, hm, af, g, t);
          }
        }
      }
    }

  }

  return;
}






#endif // LEAPFROG_H

