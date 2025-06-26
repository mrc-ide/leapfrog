#' Save parameters to HDF5 file format. This implicitly
#' processes the parameters to C++ 0 based indexing.
#'
#' @param df list/dataframe to serialize
#' @param file_path where to save the HDF5 file
#'
#' @export
save_parameters <- function(df, file_path) {
  df <- process_parameters_to_cpp(df)
  save_hdf5_file(df, file_path)
}

#' This function saves an arbitrary list/dataframe to HDF5
#' file format. Note that this file format can be used
#' in any language we support which means some R specific
#' features of the data are not guaranteed to be preserved,
#' e.g. dimnames
#'
#' @param df list/dataframe to serialize
#' @param file_path where to save the HDF5 file
#'
#' @export
save_hdf5_file <- function(df, file_path) {
  file.create(file_path)
  h5f <- hdf5r::H5File$new(file_path, mode = "w")
  invisible(tryCatch(
    save_datasets(h5f, df),
    finally = {
      h5f$close_all()
    }
  ))
}

save_datasets <- function(h5f, df, group = "") {
  lapply(names(df), function(name) {
    dat <- df[[name]]
    new_path <- sprintf("%s/%s", group, name)
    if (is.list(dat)) {
      h5f$create_group(new_path)
      save_datasets(h5f, dat, new_path)
    } else {
      h5f[[new_path]] <- dat
    }
  })
}

#' Read parameters from HDF5 file format. This implicitly
#' processes the parameters from C++ 0 based indexing to
#' R 1 based indexing.
#'
#' @param file_path HDF5 file to read
#'
#' @export
read_parameters <- function(file_path) {
  df <- read_hdf5_file(file_path)
  df <- process_parameters_to_r(df)
  df
}

#' This function reads an HDF5 file to an R dataframe.
#'
#' @param file_path where the HDF5 file is
#'
#' @return deserialized HDF5 file as a list
#'
#' @export
read_hdf5_file <- function(file_path) {
  h5f <- hdf5r::H5File$new(file_path, mode = "r+")
  groups <- hdf5r::list.groups(h5f, full.names = TRUE)
  tryCatch(
    read_datasets(h5f, groups),
    finally = {
      h5f$close_all()
    }
  )
}

read_datasets <- function(h5, groups, curr_group = "") {
  df <- lapply(names(h5), function(name) {
    new_path <- sprintf("%s/%s", curr_group, name)
    if (new_path %in% groups) {
      read_datasets(h5[[new_path]], groups, new_path)
    } else {
      obj <- h5[[new_path]]
      ret <- unlist(hdf5r::readDataSet(obj))
      dim(ret) <- obj$dims
      ret
    }
  })
  names(df) <- names(h5)
  df
}
