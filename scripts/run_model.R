#!/usr/bin/env Rscript
"Run leapfrog model and save output to specified dir
Usage:
  run_model (--output-dir=<output-dir>)

Options:
  -h --help                  Show this screen.
  --output-dir=<output-dir>  Path to save output to.
" -> usage

dat <- docopt::docopt(usage)
names(dat) <- gsub("-", "_", names(dat), fixed = TRUE)

if (!dir.exists(dat$output_dir)) {
  dir.create(dat$output_dir, recursive = TRUE)
}

input_data_dir <- frogger:::frogger_file("standalone_model", "data", "child_data")
input_data <- lapply(list.files(input_data_dir), function(filename) {
  frogger:::deserialize_tensor_to_r(file.path(input_data_dir, filename))
})
names(input_data) <- list.files(input_data_dir)

## On disk we have the C++ name saved, so we need to map to the R name
name_mapping <- frogger::frogger_input_name_mapping()
r_names <- vapply(name_mapping, "[[", character(1), "r_name")
cpp_names <- vapply(name_mapping, "[[", character(1), "cpp_name")
struct_names <- vapply(name_mapping, "[[", character(1), "struct")

name_on_disk <- sprintf("%s_%s", tolower(struct_names), cpp_names)
names(r_names) <- name_on_disk

names(input_data) <- r_names[names(input_data)]
input_data$cd4_initdist_full <- input_data[["cd4_initdist"]]
input_data$cd4_prog_full <- input_data[["cd4_prog"]]
input_data$cd4_mort_full <- input_data[["cd4_mort"]]
input_data$art_mort_full <- input_data[["art_mort"]]
input_data$t_ART_start <- 31L
input_data$projection_period <- "calendar"

out <- frogger::run_model(input_data, input_data, NULL, NULL, NULL,
                          hiv_age_stratification = "full",
                          run_child_model = TRUE)

message("Model fit complete, saving output")

for (output in names(out)) {
  frogger:::serialize_r_to_tensor(out[[output]],
                                  file.path(dat$output_dir, output))
}

message(sprintf("Saved output to dir '%s'", dat$output_dir))

