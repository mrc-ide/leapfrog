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
  int hAG, int hDS, int hDS_adol,int trans,int tx_time,int hTS>
  void leapfrog_sim(const Type *p_basepop,
                    const Type *p_sx,
                    const Type *p_netmigr,
                    const Type *p_asfr,
                    const Type *p_births_sex_prop,
                    //
                    // adult incidence
                    const Type *p_incidinput,
                    const Type *p_incrr_sex,
                    const Type *p_fert_rat,
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
                    const Type *p_scalar_art,
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
                    const double ctx_effect,
                    const Type *p_paed_art_val,
                    const int *p_artpaeds_isperc,
                    const Type *p_paed_art_elig_age,
                    const Type *p_paed_art_elig_cd4,
                
                    
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
                    Type *p_coarse_totpop1
                      ) {
    


  // macros
  // TODO: unsure if these should be defined here or elsewhere. Maybe in a future class.
  typedef Eigen::TensorMap<Eigen::Tensor<const int, 1>> TensorMapX1cI;
  typedef Eigen::TensorMap<Eigen::Tensor<const int, 2>> TensorMapX2cI;

  typedef Eigen::TensorMap<Eigen::Tensor<const Type, 1>> TensorMapX1cT;
  typedef Eigen::TensorMap<Eigen::Tensor<const Type, 2>> TensorMapX2cT;
  typedef Eigen::TensorMap<Eigen::Tensor<const Type, 3>> TensorMapX3cT;
  typedef Eigen::TensorMap<Eigen::Tensor<const Type, 4>> TensorMapX4cT;

  typedef Eigen::TensorMap<Eigen::Tensor<Type, 1>> TensorMapX1T;
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
  // mkw: hard coded number of 5 year age groups in fert rat for rn
  const TensorMapX2cT fert_rat(p_fert_rat, pAG, sim_years);
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
  const TensorMapX2cT scalar_art(p_scalar_art,  pIDX_HIVADULT, sim_years);
  const TensorMapX1cT paed_incid_input(p_paed_incid_input, sim_years);
  const TensorMapX1cT paed_cd4_dist(p_paed_cd4_dist, hDS);
  const TensorMapX1cT paed_cd4_prog(p_paed_cd4_prog, hDS);
  const TensorMapX1cT adol_cd4_prog(p_adol_cd4_prog, hDS_adol);
  const TensorMapX3cT paed_cd4_mort(p_paed_cd4_mort, hDS, trans, 5);
  const TensorMapX3cT adol_cd4_mort(p_adol_cd4_mort, hDS_adol, trans, 10);
  const TensorMapX1cT ctx_val(p_ctx_val, sim_years);
  const TensorMapX1cT paed_art_val(p_paed_art_val, sim_years);
  const TensorMapX1cT paed_art_elig_age(p_paed_art_elig_age, sim_years);
  
  const TensorMapX2cT init_art_dist(p_init_art_dist, pIDX_HIVADULT, sim_years);
  const TensorMapX3cT mort_art_rr(p_mort_art_rr, hTS, pIDX_HIVADULT, sim_years);
  

  //5 age categories
  const TensorMapX2cT paed_art_elig_cd4(p_paed_art_elig_cd4, pIDX_FERT, sim_years);
  //count by pct
  const TensorMapX2cT paed_cd4_transition(p_paed_cd4_transition, hDS_adol, hDS);

  const TensorMapX3cT paed_art_mort(p_paed_art_mort, hDS, hTS, 5);
  const TensorMapX3cT adol_art_mort(p_adol_art_mort, hDS_adol, hTS, 10);
  const TensorMapX1cI artpaeds_isperc(p_artpaeds_isperc, sim_years);
  


  // outputs
  TensorMapX3T totpop1(p_totpop1, pAG, NG, sim_years);
  TensorMapX3T hivpop1(p_hivpop1, pAG, NG, sim_years);
  TensorMapX3T hivnpop1(p_hivnpop1, pAG, NG, sim_years);
  TensorMapX3T infections(p_infections, pAG, NG, sim_years);
  TensorMapX1T births(p_births, sim_years); 
  TensorMapX1T hiv_births(p_hiv_births, sim_years); 
  TensorMapX3T natdeaths(p_natdeaths, pAG, NG, sim_years);
  TensorMapX3T natdeaths_hivpop(p_natdeaths_hivpop, pAG, NG, sim_years);
  TensorMapX3T hivdeaths(p_hivdeaths, pAG, NG, sim_years);
  TensorMapX4T aidsdeaths_noart(p_aidsdeaths_noart, hDS, hAG, NG, sim_years);
  TensorMapX5T aidsdeaths_art(p_aidsdeaths_art, hTS, hDS, hAG, NG, sim_years);
  TensorMapX4T artinit(p_artinit, hDS, hAG, NG, sim_years);

  TensorMapX4T hivstrat_adult(p_hivstrat_adult, hDS, hAG, NG, sim_years);
  TensorMapX5T artstrat_adult(p_artstrat_adult, hTS, hDS, hAG, NG, sim_years);
  //adding in transmission category
  TensorMapX5T hivstrat_paeds(p_hivstrat_paeds, hDS, trans, pIDX_HIVADULT, NG, sim_years);
  TensorMapX5T artstrat_paeds(p_artstrat_paeds, hTS, hDS, pIDX_HIVADULT, NG, sim_years);
  TensorMapX5T artelig_paeds(p_artelig_paeds, hDS, trans, pIDX_HIVADULT, NG, sim_years);
  TensorMapX4T aidsdeaths_noart_paed(p_aidsdeaths_noart_paed, hDS, pIDX_HIVADULT, NG, sim_years);
  TensorMapX5T aidsdeaths_art_paed(p_aidsdeaths_art_paed, hTS, hDS, pIDX_HIVADULT, NG, sim_years);
  TensorMapX5T deaths_paeds(p_deaths_paeds, hDS, trans, pIDX_HIVADULT, NG, sim_years);  
  TensorMapX5T deaths_paeds_art(p_deaths_paeds_art, hDS, trans, pIDX_HIVADULT, NG, sim_years);  
  TensorMapX5T grad_paeds(p_grad_paeds, hDS, trans, pIDX_HIVADULT, NG, sim_years);  
  TensorMapX5T grad_paeds_art(p_grad_paeds_art, hTS, hDS,  pIDX_HIVADULT, NG, sim_years);  
  TensorMapX1T artnum_paed(p_artnum_paed, sim_years); 
  TensorMapX5T init_art_paed(p_init_art_paed, hDS, trans, pIDX_HIVADULT, NG, sim_years);  
  

  TensorMapX3T coarse_totpop1(p_coarse_totpop1, hAG, NG, sim_years);
  
  
  
  

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
  
  
  
  int everARTelig_idx = hDS;

  ////////////////////////////////////
  ////  do population projection  ////
  ////////////////////////////////////

  for(int t = 1; t < sim_years; t++){

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
      int a = pIDX_HIVADULT;
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
          if(t > t_ART_start)
            for(int hu = 0; hu < hTS; hu++) {
              artstrat_adult(hu, hm, ha, g, t) = (1.0 - hiv_ag_prob(ha, g)) * artstrat_adult(hu, hm, ha, g, t-1);
              artstrat_adult(hu, hm, ha, g, t) += hiv_ag_prob(ha-1, g) * artstrat_adult(hu, hm, ha-1, g, t-1);
            }
        }
      }
    }
    
    

    
    for(int g = 0; g < NG; g++){
      for(int ha = 1; ha < 5; ha++) {
        for(int hm = 0; hm < hDS; hm++){
          for(int cat = 0 ; cat < 4; cat++){
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
          for(int cat = 0 ; cat < trans; cat++){
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
          for(int cat = 0 ; cat < trans; cat++){
           hivstrat_paeds(hm, cat, ha, g, t) += hivstrat_paeds(hm, cat, ha-1, g, t-1) * sx(ha, g, t);
            
          }
          for(int dur = 0; dur < hTS; dur++){
            artstrat_paeds(dur, hm, ha, g, t) += artstrat_paeds(dur, hm, ha-1, g, t-1) * sx(ha, g, t);
          }
            
        }
      }
    }
    
    //Progress ART to the correct time on ART
    for(int hm = 0; hm < hDS; hm++){
      for(int af = 0; af < pIDX_FERT; af++){
        for(int g = 0; g < NG; g++){
       //   artstrat_paeds(2, hm, af, g, t) += artstrat_paeds(0, hm, af, g, t) > 0 ? artstrat_paeds(0, hm, af, g, t)  : 0;
      //    artstrat_paeds(0, hm, af, g, t) -= artstrat_paeds(0, hm, af, g, t) > 0 ? artstrat_paeds(0, hm, af, g, t) : 0;
        }
      }
    }
    
    // !!!TODO: add HIV+ 15 year old entrants
    for (int g = 0; g < NG; g++) {
      for (int hm = 0; hm < hDS; hm++) {
        hivstrat_adult(hm, 0, g, t) = (1.0 - hiv_ag_prob(0, g)) * hivstrat_adult(hm, 0, g, t-1);
        hivstrat_adult(hm, 0, g, t) = hivstrat_adult(hm, 0, g, t) + age15hivpop(1, hm, g, t);
        // ADD HIV+ entrants here
        if(t > t_ART_start) {
          for(int hu = 0; hu < hTS; hu++) {
            artstrat_adult(hu, hm, 0, g, t) = (1.0 - hiv_ag_prob(0, g)) * artstrat_adult(hu, hm, 0, g, t-1);
            // ADD HIV+ entrants here
            //       artpop_t(hu, hm, 0, g, t) += paedsurv_g * paedsurv_artcd4dist(hu, hm, g, t) * entrantartcov(g, t);
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
          if(t > t_ART_start) {
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
    // loop hiv_steps_per_year
    // adjust population to match target population size
    // TensorFixedSize<Type, Sizes<hAG, NG>> popadjprob;
    // TensorFixedSize<Type, Sizes<hAG, NG>> hivpopadjprob;
    for(int g = 0; g < NG; g++){
      for(int ha = 0; ha < pAG; ha++){
        //popadjprob(ha, g) = 0.0;
        //popadjprob(ha, g) = basepop(ha, g, t) / totpop1(ha, g, t);
        // maybe need to change this for coarse age groups
        // hivpopadjprob(ha, g) = popadjprob(ha, g) ;
        
        totpop1(ha, g, t) = basepop(ha, g, t);
        hivpop1(ha, g, t) =  hivpop1(ha, g, t) * basepop(ha, g, t) / totpop1(ha, g, t);
        
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
       hivnpop1(a, g, t) = (basepop(a, g, t)) - hivpop1(a, g, t);
        
      }
    }
    hiv_births(t) = 0.0;
    for(int af = 0; af < pAG; af++) {
      //hiv_births(t) += (hivpop1(pIDX_FERT + af, FEMALE, t + hivpop1(pIDX_FERT + af, FEMALE, t-1) * 0.5) * asfr(af, t) * fert_rat(pIDX_FERT + af, t));
    }



    for(int g = 0; g < NG; g++){
      //WAS 6
      for(int hm = 0; hm < hDS; hm++){
        for(int af = 1; af < 5; af++){
          for(int cat = 0; cat < trans; cat++){
            deaths_paeds(hm, cat, af, g, t) += hivstrat_paeds(hm, cat, af, g, t) - (1 - ctx_effect * ctx_val(t)) * hivstrat_paeds(hm, cat, af, g, t) * paed_cd4_mort(hm, cat, af); 
            aidsdeaths_noart_paed(hm, af, g, t) +=  (1 - ctx_effect * ctx_val(t)) * hivstrat_paeds(hm, cat, af, g, t) * paed_cd4_mort(hm, cat, af); // output hiv deaths, aggregated across transmission category
          }
        }
      }
    }
    
    for(int g = 0; g < NG; g++){
      for(int hm = 0; hm < hDS_adol; hm++){
        for(int af = 5; af < pIDX_HIVADULT; af++){
          for(int cat = 0; cat < trans; cat++){
            //deaths_paeds kinda a misnomer, those alive after hiv mortality
             deaths_paeds(hm, cat, af, g, t) += hivstrat_paeds(hm, cat, af, g, t) ;
             deaths_paeds(hm, cat, af, g, t) -= (1 - ctx_effect * ctx_val(t)) * hivstrat_paeds(hm, cat, af, g, t) * adol_cd4_mort(hm, cat, af - 5); 
             aidsdeaths_noart_paed(hm, af, g, t) +=  (1 - ctx_effect * ctx_val(t)) * hivstrat_paeds(hm, cat, af, g, t) * adol_cd4_mort(hm, cat, af - 5); // output hiv deaths, aggregated across transmission category
             
          }
        }
      }
    }
    
    

    //progress through CD4 categories
    for(int g = 0; g < NG; g++){
      //the nosocomial infections aren't distributed so can't just move everything forward. So going to limit hm to just go to the n+1 basically
      //WAS 6
      for(int hm = 1; hm < hDS; hm++){
        for(int af = 1; af < 5; af++){
          for(int cat = 0; cat < trans; cat++){
            grad_paeds(hm - 1, cat, af, g, t) -= ((hivstrat_paeds(hm-1, cat, af-1, g, t-1) + artstrat_paeds(hm-1, cat, af-1, g, t-1)) > 0) ? (deaths_paeds(hm - 1, cat, af, g, t) * paed_cd4_prog(hm - 1) + hivstrat_paeds(hm - 1, cat, af, g, t) * paed_cd4_prog(hm - 1)) / 2: 0.0; //moving to next cd4 category
            grad_paeds(hm, cat, af, g, t) += ((hivstrat_paeds(hm-1, cat, af-1, g, t-1) + artstrat_paeds(hm-1, cat, af-1, g, t-1)) > 0) ? (deaths_paeds(hm - 1, cat, af, g, t) * paed_cd4_prog(hm - 1) + hivstrat_paeds(hm - 1, cat, af, g, t) * paed_cd4_prog(hm - 1)) / 2 : 0.0; //moving into this cd4 category
          }
        }
      }
    }
    
    //progress through CD4 categories
    //changed here
    for(int g = 0; g < NG; g++){
      for(int hm = 1; hm < hDS_adol; hm++){
        for(int af = 5; af < pIDX_FERT; af++){
          for(int cat = 0; cat < trans; cat++){
            //deaths_paeds causing NANs
            grad_paeds(hm - 1, cat, af, g, t) -= (deaths_paeds(hm - 1, cat, af, g, t) * adol_cd4_prog(hm - 1) + hivstrat_paeds(hm - 1, cat, af, g, t) * adol_cd4_prog(hm - 1)) / 2; //moving to next cd4 category
            grad_paeds(hm, cat, af, g, t) += (deaths_paeds(hm - 1, cat, af, g, t) * adol_cd4_prog(hm - 1) + hivstrat_paeds(hm - 1, cat, af, g, t) * adol_cd4_prog(hm - 1)) / 2; //moving into this cd4 category
          }
        }
      }
    }

    
    for(int g = 0; g < NG; g++){
      //WAS 6
      for(int hm = 0; hm < hDS; hm++){
        for(int af = 1; af < 5; af++){
          for(int cat = 0; cat < trans; cat++){
            grad_paeds(hm, cat, af, g, t) -= hivstrat_paeds(hm, cat, af, g, t) * paed_cd4_mort(hm, cat, af) * ctx_val(t) * (1 - ctx_effect) + hivstrat_paeds(hm, cat, af, g, t) * paed_cd4_mort(hm, cat, af) * (1 - ctx_val(t));
            hivstrat_paeds(hm, cat, af, g, t) += grad_paeds(hm, cat, af, g, t) ; 
          }

        }
      }
    }
    
    //not accounting for ctx 5 year 
    for(int g = 0; g < NG; g++){
      for(int hm = 0; hm < hDS_adol; hm++){
        for(int af = 5; af < pIDX_FERT; af++){
          for(int cat = 0; cat < trans; cat++){
           grad_paeds(hm, cat, af, g, t) -= hivstrat_paeds(hm, cat, af, g, t) * adol_cd4_mort(hm, cat, af - 5) * ctx_val(t) * (1 - ctx_effect) + hivstrat_paeds(hm, cat, af, g, t) * adol_cd4_mort(hm, cat, af - 5) * (1 - ctx_val(t));
           hivstrat_paeds(hm, cat, af, g, t) += grad_paeds(hm, cat, af, g, t) ; 
          }
        }
      }
    }
    
    
    TensorFixedSize<Type, Sizes<hDS, trans, pIDX_HIVADULT, NG>> need_art_paed;
    need_art_paed.setZero();

 
    
    //all children under a certain age eligible for ART
    for(int g = 0; g < NG; g++){
      for(int cat = 0; cat < trans; cat++){
        for(int af = 0; af < pIDX_HIVADULT; af++){
          for(int hm = 0; hm < hDS; hm++){
            need_art_paed(hm, cat, af, g) += af < paed_art_elig_age(t) ? hivstrat_paeds(hm, cat, af, g, t) :0.0 ;
          }
        } 
      }
    }
    
    
    //all children under a certain CD4 are eligible for ART, regardless of age
    for(int g = 0; g < NG; g++){
      for(int cat = 0; cat < trans; cat++){
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
            for(int cat = 0; cat < trans; cat++){
              artnum_paed(t) += need_art_paed(hm, cat, af, g);
            }
           for(int dur = 0; dur < hTS; dur++){
            //  artnum_paed(t) += artstrat_paeds(dur, hm, af, g, t);
          }
        } 
      }
    }
    
   // if ( (!art15plus_isperc(g, t-2)) & (!art15plus_isperc(g, t-1)) ){ // both numbers
  //    artnum_hts = (0.5-dt*(hts+1))*art15plus_num(g, t-2) + (dt*(hts+1)+0.5)*art15plus_num(g, t-1);
  //  } else if (art15plus_isperc(g, t-2) & art15plus_isperc(g, t-1)){ // both percentages
   //   Type artcov_hts = (0.5-dt*(hts+1))*art15plus_num(g, t-2) + (dt*(hts+1)+0.5)*art15plus_num(g, t-1);
   //   artnum_hts = artcov_hts * (Xart_15plus + Xartelig_15plus);
   // } else if ( (!art15plus_isperc(g, t-2)) & art15plus_isperc(g, t-1)) { // transition from number to percentage
  //    Type curr_coverage = Xart_15plus / (Xart_15plus + Xartelig_15plus);
 //     Type artcov_hts = curr_coverage + (art15plus_num(g, t-1) - curr_coverage) * dt / (0.5-dt*hts);
 //     artnum_hts = artcov_hts * (Xart_15plus + Xartelig_15plus);
//    }



    
    //how many should initialize ART
    for(int g = 0; g < NG; g++){
        for(int af = 0; af < pIDX_HIVADULT; af++){
          for(int hm = 0; hm < hDS; hm++){
            for(int cat = 0; cat < trans; cat++){
              init_art_paed(hm, cat, af, g, t) += need_art_paed(hm, cat, af, g);
              for(int dur = 0; dur < hTS; dur++){
            //    init_art_paed(hm, cat, af, g, t) -= artstrat_paeds(dur, hm, af, g, t) / 4;
                init_art_paed(hm, cat, af, g, t) = init_art_paed(hm, cat, af, g, t) < 0 ? 0.0 : init_art_paed(hm, cat, af, g, t);
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
            for(int cat = 0; cat < trans; cat++){
            init_art_paed_total += init_art_paed(hm, cat, af, g, t);
          }
        } 
      }
    }
    
    
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

    
    //Progress ART to the correct time on ART
    for(int hm = 0; hm < hDS; hm++){
      for(int af = 0; af < pIDX_FERT; af++){
        for(int g = 0; g < NG; g++){
          artstrat_paeds(1, hm, af, g, t) += artstrat_paeds(0, hm, af, g, t) > 0 ? artstrat_paeds(0, hm, af, g, t) : 0;
          artstrat_paeds(0, hm, af, g, t) -= artstrat_paeds(0, hm, af, g, t) > 0 ? artstrat_paeds(0, hm, af, g, t) : 0;
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
    
    if ( (!artpaeds_isperc(t)) & (!artpaeds_isperc(t-1)) ){ // both numbers
      artnum_paed(t) = paed_art_val(t-1) > 0 ?  (paed_art_val(t) + paed_art_val(t-1)) / 2 : 0.0;
      //Remove how many that are already on ART
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
      artnum_paed(t) = paed_art_val(t-1) > 0 ? artnum_paed(t) * (paed_art_val(t) + paed_art_val(t-1)) / 2 : 0.0;
      //Remove how many that are already on ART
      for(int g = 0; g < NG; g++){
        for(int af = 0; af < pIDX_HIVADULT; af++){
          for(int hm = 0; hm < hDS; hm++){
            for(int dur = 0; dur < hTS; dur++){
             // artnum_paed(t) -= (artstrat_paeds(dur, hm, af, g, t) )  ;
              
            }
          }
        }
      }
      artnum_paed(t) = init_art_paed_total < artnum_paed(t) ? init_art_paed_total : artnum_paed(t) ;
      //if ART coverage in percentages isn't changing we don't continue to add more people on art (idrk get this)
      artnum_paed(t) = (paed_art_val(t-1) >= paed_art_val(t)) & paed_art_val(t) != 1 ? 0 : artnum_paed(t); 
      


    } else if (artpaeds_isperc(t) & !artpaeds_isperc(t-1)){ // num to percentage
      
      //Remove how many that are already on ART
      double temp ;
      temp = 0.0;
      for(int g = 0; g < NG; g++){
        for(int af = 0; af < pIDX_HIVADULT; af++){
          for(int hm = 0; hm < hDS; hm++){
            for(int dur = 0; dur < hTS; dur++){
              artnum_paed(t) += artstrat_paeds(dur, hm, af, g, t)   ;
              
            }
          }
        }
      }

      artnum_paed(t) = (paed_art_val(t-1) + (artnum_paed(t)) * paed_art_val(t)) / 2 ;

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
      artnum_paed(t) = init_art_paed_total < paed_art_val(t) ? init_art_paed_total : paed_art_val(t) ;
      
      double last_year_art;
      last_year_art = 0.0;
      for(int g = 0; g < NG; g++){
        for(int af = 0; af < pIDX_HIVADULT; af++){
          for(int hm = 0; hm < hDS; hm++){
            for(int dur = 0; dur < hTS; dur++){
              last_year_art += (artstrat_paeds(dur, hm, af, g, t) )  ;
              
            }
          }
        }
      }
   
      artnum_paed(t) = (last_year_art + artnum_paed(t)) / 2 ;
      
      //Remove how many that are already on ART
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
      
      
    }
    
    
    double initByAge = 0.0;
    for(int g = 0; g < NG; g++){
      for(int af = 0; af < pIDX_HIVADULT; af++){
        for(int hm = 0; hm < hDS; hm++){
          for(int cat = 0; cat < trans; cat++){
            initByAge +=  init_art_paed(hm, cat, af, g, t) * init_art_dist(af, t);
          }
        }
      }
    }
    double adj ;
    adj = initByAge == 0 ? 1 : artnum_paed(t) / initByAge ;
    
    //the nosocomial infections aren't distributed so can't just move everything forward. So going to limit hm to just go to the n+1 basically
    for(int g = 0; g < NG; g++){
      for(int cat = 0; cat < trans; cat++){
        for(int af = 0; af < pIDX_HIVADULT; af++){
          for(int hm = 0; hm < hDS; hm++){
            double scalar ;
            scalar = (adj * init_art_dist(af, t)) > 1 ? 1 : adj * init_art_dist(af, t);
            scalar = artnum_paed(t) > 0 ? scalar : 0.0;
            //not really sure why this needs to be here, but without it it scales down certain age groups when there is enough ART
            scalar = (paed_art_mort(hm, cat, af) == 0) ?  1 : scalar;
          
            artstrat_paeds(0, hm, af, g, t) +=  scalar * init_art_paed(hm, cat, af, g, t) ;
         
            hivstrat_paeds(hm, cat, af, g, t) -=  scalar* init_art_paed(hm, cat, af, g, t) ;
            
          }
        }
      }
    }
    
    
    
    for(int g = 0; g < NG; g++){
      for(int hm = 0; hm < hDS; hm++){
        for(int af = 0; af < pIDX_HIVADULT; af++){
          //SOMETHING HAPPENING HERE
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

  
    
    //distribute across eligible ages, right now just going to hardcode
    for(int g = 0; g < NG; g++){
      for(int af = 0; af < 5; af++){
        infections(af, g, t) = 0 ;
        infections(af, g, t) = paed_incid_input(t) / 10;
        hivpop1(af, g, t) += infections(af, g, t);
 
        for(int hm = 0; hm < hDS; hm++){
          for(int cat = 0; cat < 1; cat++){
            //putting them all in perinatal trans to match spec nosocomial
            hivstrat_paeds(hm, cat, af, g, t) += paed_cd4_dist(hm) > 0 ? infections(af, g, t) * paed_cd4_dist(hm) : 0.0;
          }
        }
      }
    }
    
    
  }

  return;
}






#endif // LEAPFROG_H
