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
  pmtct_dropout[is.na(pmtct_dropout)] <- 0
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
pmtct_input_isperc <- array(data = 0, dim = c(nrow(pmtct[[1]]),ncol(pmtct[[1]])), dimnames = list(hp = unlist(lapply(rownames(pmtct[[1]]), gsub, pattern = '_number', replacement = '')),
                                                                                                  year = 1970:(1969 + ncol(pmtct[[1]]))))
pmtct_input_isperc[!is.na(pmtct[[1]])] <- 1

pmtct_input_isperc
}

input_pmtct_retained <- function(pjnz){
  pmtct_retained <- dp_read_pmtct_retained(pjnz)
  dropout <- format_pmtct(pjnz)$dropout

  ##ASK: not sure of the best way to do these indexes
  out <- array(data = NA, dim = c(6, ncol(pmtct_retained)), dimnames = list(type = c(rownames(pmtct_retained), rownames(dropout)), year = colnames(pmtct_retained)))
  out[1:2,] <- pmtct_retained / 100
  out[3:nrow(out),] <- unlist(dropout)

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

input_mothers_reallocated <- function(pjnz){
  patients_reallocated <- dp_read_mothers_reallocated(pjnz)
  out <- array(patients_reallocated, dim = length(patients_reallocated), dimnames = list(year = names(patients_reallocated)))

  out
}

input_breastfeeding_dur <- function(pjnz){
  bf <- dp_read_breastfeeding(pjnz)
  dimnames(bf)

  dpfile <- grep(".DP$", utils::unzip(pjnz, list=TRUE)$Name, value=TRUE)
  dp <- utils::read.csv(unz(pjnz, dpfile), as.is=TRUE)
  dpsub_new <- function(tag, rows, cols, tagcol=1){
    dp[which(dp[,tagcol]==tag)+rows, cols]
  }

  yr_start <- as.integer(dpsub_new("<FirstYear MV2>",2,4))
  yr_end <- as.integer(dpsub_new("<FinalYear MV2>",2,4))
  proj.years <- yr_start:yr_end
  timedat.idx <- 4+1:length(proj.years)-1

  out <- array(data = NA, dim = c(nrow(bf$notbreastfeeding_percent_noarv), ncol(bf$notbreastfeeding_percent_noarv), 2),
               dimnames = list(child_age_months = rownames(bf$notbreastfeeding_percent_arv), year = proj.years, trt = c('no_art', 'art')))
  out[,,1] <- bf$notbreastfeeding_percent_noarv / 100
  out[,,2] <- bf$notbreastfeeding_percent_arv / 100

  out

}

input_childart <- function(pjnz){
  art <- dp_read_childart(pjnz)

  cotrim <- array(art$childart[1,], dim = c(1, ncol(art$childart)), dimnames = list(val = '', year = colnames(art$childart)))
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

















