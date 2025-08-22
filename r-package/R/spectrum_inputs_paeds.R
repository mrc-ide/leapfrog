###############################################
###Paediatric input formatting
###############################################
format_pmtct <- function(pjnz){
  pmtct_arv <- data.frame(dp_read_pmtct(pjnz))

  ## remove all totals
  pmtct_arv_labels <- c("pmtct_noprophylaxis_percent",
                        "pmtct_total_number",
                        "pmtct_postnatal_noprophylaxis_percent",
                        "pmtct_postnatal_total_number",
                        "pmtct_postnatal_total_number")
  pmtct_arv_rows_keep <- !(rownames(pmtct_arv) %in% pmtct_arv_labels)
  pmtct_arv <- pmtct_arv[pmtct_arv_rows_keep, ]

  ##separate out dropouts
  pmtct_dropout <- pmtct_arv %>% dplyr::filter(dplyr::row_number() %in% grep(pattern = 'dropout', rownames(pmtct_arv))) / 100
  pmtct_dropout[is.na(pmtct_dropout)] <- 0
  pmtct_arv <- pmtct_arv %>% dplyr::filter(!dplyr::row_number() %in% grep(pattern = 'dropout', rownames(pmtct_arv)))

  ## take out postnatal
  pmtct_arv_postnatal <- pmtct_arv %>% dplyr::filter(dplyr::row_number() %in% grep(pattern = 'postnatal', rownames(pmtct_arv)))
  pmtct_arv <- pmtct_arv %>% dplyr::filter(!dplyr::row_number() %in% grep(pattern = 'postnatal', rownames(pmtct_arv)))

  ## number v percents
  pmtct_arv_number <- pmtct_arv %>% dplyr::filter(dplyr::row_number() %in% grep(pattern = 'number', rownames(pmtct_arv)))
  pmtct_arv_percent <- pmtct_arv %>% dplyr::filter(dplyr::row_number() %in% grep(pattern = 'percent', rownames(pmtct_arv)))


  if(all(is.na(pmtct_arv_number))){
    pmtct_number <- F
    pmtct_val = pmtct_arv_percent / 100
  }else{
    pmtct_number <- T
    pmtct_val = pmtct_arv_number
  }

  list(pmtct_value = list(num = pmtct_arv_number, pct = pmtct_arv_percent),
       pmtct_number = pmtct_number,
       dropout = pmtct_dropout,
       postnatal = pmtct_arv_postnatal)
}

input_pmtct <- function(pjnz){
  pmtct <- format_pmtct(pjnz)$pmtct_value
  out <- array(data = NA, dim = c(nrow(pmtct[[1]]), ncol(pmtct[[1]]), 2), dimnames = list(type = rownames(pmtct[[1]]), year = 1970:(1969 + ncol(pmtct[[1]])), metric = c('number', 'percent')))
  out[,,1] <- unlist(pmtct$num)
  out[,,2] <- unlist(pmtct$pct)
  out[is.na(out)] <- 0

  out
}

input_pmtct_ispercent <- function(pjnz){

pmtct <- format_pmtct(pjnz)$pmtct_value

## pmtct_input_isperc
pmtct_input_isperc <- array(data = FALSE, dim = c(nrow(pmtct[[1]]),ncol(pmtct[[1]])), dimnames = list(hp = unlist(lapply(rownames(pmtct[[1]]), gsub, pattern = '_number', replacement = '')),
                                                                                                  year = 1970:(1969 + ncol(pmtct[[1]]))))
pmtct_input_isperc[!is.na(pmtct[[1]])] <- TRUE

pmtct_input_isperc
}

input_pmtct_retained <- function(pjnz){
  pmtct_retained <- dp_read_pmtct_retained(pjnz)
  dropout <- format_pmtct(pjnz)$dropout

  ##ASK: not sure of the best way to do these indexes
  out <- array(data = 0, dim = c(7, 3, ncol(pmtct_retained)),
               dimnames = list(pmtct = c('Option A', 'Option B', 'SDNVP',
                                         'Dual ARV', 'Option B+: before pregnancy',
                                         'Option B+: >4 weeks', 'Option B+: <4 weeks'),
                               trans_time = c('perinatal',
                                              'breastfeeding <12 months',
                                              'breastfeeding 12+ months'),
                               year = colnames(pmtct_retained)))
  out[,'perinatal',] <- 1
  out[c('Option B+: before pregnancy', 'Option B+: >4 weeks'),'perinatal',] <- pmtct_retained / 100
  out[c('Option A', 'Option B','Option B+: before pregnancy', 'Option B+: >4 weeks','Option B+: <4 weeks'),'breastfeeding <12 months',] <- matrix(unlist(dropout['pmtct_postnatal_monthlydropout_art0to12months',]), nrow = 5, ncol = 61, byrow = TRUE)
  out[c('Option A', 'Option B','Option B+: before pregnancy', 'Option B+: >4 weeks','Option B+: <4 weeks'),'breastfeeding 12+ months',] <- matrix(unlist(dropout['pmtct_postnatal_monthlydropout_art12plusmonths',]), nrow = 5, ncol = 61, byrow = TRUE)

  out
}

input_abortion <- function(pjnz){
  abort <- dp_read_abortion(pjnz)
  abort$pregtermabortion_ispercent[is.na(abort$pregtermabortion_ispercent)] <- TRUE
  abort$pregtermabortion[abort$pregtermabortion_ispercent] <- abort$pregtermabortion[abort$pregtermabortion_ispercent] / 100

  out <- array(data = NA, dim = c(2, length(abort$pregtermabortion)), dimnames = list(type = c('value', 'ispercent'), year = names(abort$pregtermabortion)))
  out[1,] <- abort$pregtermabortion
  ##this is changing to 0/1
  out[2,] <- abort$pregtermabortion_ispercent

  out
}

input_mothers_reallocated <- function(pjnz){
  patients_reallocated <- dp_read_mothers_reallocated(pjnz)
  out <- array(patients_reallocated, dim = length(patients_reallocated), dimnames = list(year = names(patients_reallocated)))

  out
}

input_breastfeeding_dur <- function(pjnz){
  bf <- dp_read_breastfeeding(pjnz)
  dimnames(bf)

  out <- array(data = NA, dim = c(nrow(bf$notbreastfeeding_percent_noarv), ncol(bf$notbreastfeeding_percent_noarv), 2), dimnames = list(child_age_months = rownames(bf$notbreastfeeding_percent_arv), year = (1970:2030)[1:ncol(bf$notbreastfeeding_percent_noarv)], trt = c('no_art', 'art')))
  out[,,1] <- bf$notbreastfeeding_percent_noarv / 100
  out[,,2] <- bf$notbreastfeeding_percent_arv / 100

  out

}

input_childart <- function(pjnz){
  art <- dp_read_childart(pjnz)

  cotrim <- array(art$childart[1,], dim = c(ncol(art$childart)), dimnames = list(year = colnames(art$childart)))
  cotrim_pct <- array(art$childart_ispercent[1,], dim = c(1, ncol(art$childart)), dimnames = list(val = '', year = colnames(art$childart)))
  ######################
  ##Sort whether ART is being input by age or aggregate
  ######################
  age_specific <- rep(F, ncol(art$childart))
  if(any(which(!is.na(art$childart[3,])))){
    age_specific[which(!is.na(art$childart[3,]))] = T
  }else{
    age_specific[] = F
  }

  out <- art$childart[2:5,]
  out[is.na(out)] <- 0

  list(ctx = cotrim,
       ctx_percent = cotrim_pct,
       age_spec = age_specific,
       child_art = out,
       art_ispercent = art$childart_ispercent)

}


input_childart_ltfu <- function(pjnz){
  ltfu <- dp_read_childltfu(pjnz)
  out <- array(ltfu, dim = length(ltfu), dimnames = list(year = names(ltfu)))

  out
}

















