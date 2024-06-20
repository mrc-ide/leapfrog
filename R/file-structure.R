get_coarse_ages_lfrog <- function(input){
  sample_fp <- eppasm::prepare_directincid(test_path("../testdata/spectrum/v6.13/bwa_demproj-only_spectrum-v6.13_2022-02-12.PJNZ"))
  age_mapping <- data.frame(
    eppasm_age = sample_fp$ss$ag.idx,
    leapfrog_age = 1:66
  )
  input = as.list(plyr::alply(input, 3))
  as.data.frame.table(input, responseName = "EPPASM")
}
