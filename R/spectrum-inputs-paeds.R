###############################################
###Paediatric input formatting
###############################################
library(dplyr)

format_pmtct <- function(pjnz){
  pmtct_arv <- dp_read_pmtct(pjnz)
  
  ## remove all totals
  pmtct_arv <- pmtct_arv %>% data.frame() %>% filter(!row_number() %in% which(rownames(pmtct_arv) %in% c('pmtct_noprophylaxis_percent', 'pmtct_total_number',
                                                                                                         'pmtct_postnatal_noprophylaxis_percent', 'pmtct_postnatal_total_number',
                                                                                                         'pmtct_postnatal_total_number')))
  
  ##separate out dropouts
  pmtct_dropout <- pmtct_arv %>% data.frame() %>% filter(row_number() %in% grep(pattern = 'dropout', rownames(pmtct_arv))) / 100
  pmtct_arv <- pmtct_arv %>% data.frame() %>% filter(!row_number() %in% grep(pattern = 'dropout', rownames(pmtct_arv)))
  
  ## take out postnatal
  pmtct_arv_postnatal <- pmtct_arv %>% data.frame() %>% filter(row_number() %in% grep(pattern = 'postnatal', rownames(pmtct_arv)))
  pmtct_arv<- pmtct_arv %>% data.frame() %>% filter(!row_number() %in% grep(pattern = 'postnatal', rownames(pmtct_arv)))
  
  
  ## number v percents
  pmtct_arv_number <- pmtct_arv %>% data.frame() %>% filter(row_number() %in% grep(pattern = 'number', rownames(pmtct_arv)))
  pmtct_arv_percent <- pmtct_arv %>% data.frame() %>% filter(row_number() %in% grep(pattern = 'percent', rownames(pmtct_arv)))
  
  if(all(is.na(pmtct_arv_number))){
    pmtct_number <- F
    pmtct_val = pmtct_arv_percent / 100
  }else{
    pmtct_number <- T
    pmtct_val = pmtct_arv_number
  }
  
  list(pmtct_value = pmtct_val, 
       pmtct_number = pmtct_number,
       dropout = pmtct_dropout, 
       postnatal = pmtct_arv_postnatal)
}

input_pmtct <- function(pjnz){
  pmtct <- format_pmtct(pjnz)$pmtct_value
  out <- array(data = NA, dim = c(nrow(pmtct), ncol(pmtct)), dimnames = list(type = rownames(pmtct), year = 1970:2030))
  out[] <- unlist(pmtct)
  out[is.na(out)] <- 0
  
  out
}

input_pmtct_retained <- function(pjnz){
  pmtct_retained <- dp_read_pmtct_retained(pjnz)
  dropout <- format_pmtct(pjnz)$dropout
  
  ##ASK: not sure of the best way to do these indexes
  out <- array(data = NA, dim = c(6, ncol(pmtct_retained)), dimnames = list(type = c(rownames(pmtct_retained), rownames(dropout)), year = colnames(pmtct_retained)))
  out[1:2,] <- pmtct_retained / 100
  out[3:nrow(out),] <- unlist(dropout) / 100
  
  out
}

input_abortion <- function(pjnz){
  abort <- dp_read_abortion(pjnz)
  abort$pregtermabortion_ispercent[is.na(abort$pregtermabortion_ispercent)] <- TRUE
  abort$pregtermabortion[abort$pregtermabortion_ispercent] <- abort$pregtermabortion[abort$pregtermabortion_ispercent] / 100
  
  out <- array(data = NA, dim = c(1, length(abort$pregtermabortion), 2), dimnames = list(val = '', year = names(abort$pregtermabortion), type = c('value', 'ispercent')))
  out[,,1] <- abort$pregtermabortion
  ##this is changing to 0/1 
  out[,,2] <- abort$pregtermabortion_ispercent
  
  out
}

input_breastfeeding_dur <- function(pjnz){
  bf <- dp_read_breastfeeding(pjnz)
  dimnames(bf)
  
  out <- array(data = NA, dim = c(nrow(bf$notbreastfeeding_percent_noarv), ncol(bf$notbreastfeeding_percent_noarv), 2), dimnames = list(child_age_months = rownames(bf$notbreastfeeding_percent_arv), year = 1970:2030, trt = c('no_art', 'art')))
  out[,,1] <- bf$notbreastfeeding_percent_noarv / 100
  out[,,2] <- bf$notbreastfeeding_percent_arv / 100
  
  out
  
}

input_childart <- function(pjnz){
  art <- dp_read_childart(pjnz)
  
  cotrim <- array(art$childart[1,], dim = c(1, ncol(art$childart)), dimnames = list(val = '', year = colnames(art$childart)))
  
  ######################
  ##Sort whether ART is being input by age or aggregate
  ######################
  if(all(is.na(art$childart[3:5,]))){
    age_specific = F
  }else{
    age_specific = T
  }
  
  if(age_specific){
    out <- array(data = NA, dim = c(3,ncol(art$childart)), dimnames = list(age = row.names(art$childart)[3:5], year = colnames(art$childart)))
    out[] <- art$childart[3:5,]
  }else{
    out <- array(data = NA, dim = c(1,ncol(art$childart)), dimnames = list(age = row.names(art$childart)[2], year = colnames(art$childart)))
    out[] <- art$childart[2,]
  }
  
  list(ctx = cotrim,
       age_spec = age_specific,
       child_art = out)
  
}


input_childart_ltfu <- function(pjnz){
  ltfu <- dp_read_childltfu(pjnz)
  out <- array(ltfu, dim = length(ltfu), dimnames = list(year = names(ltfu)))
  
  out
}

















