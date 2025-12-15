get_dim_vars <- function(dp) {
  years_cfg <- get_years_cfg()
  start <- get_data_from_cfg("first_year", years_cfg$first_year, dim_vars, dp)$data
  end <- get_data_from_cfg("final_year", years_cfg$final_year, dim_vars, dp)$data

  dim_vars <- get_static_dim_vars()
  dim_vars$years <- as.character(start:end)
  dim_vars
}

get_data_from_cfg <- function(name, cfg, dim_vars, dp) {
  data <- NULL
  matching_tag_cfg <- NULL
  for (tag_cfg in cfg$read) {
    data <- get_data_from_tag_cfg(tag_cfg, dim_vars, dp)
    matching_tag_cfg <- tag_cfg
    if (!is.null(data)) break()
  }
  if (is.null(data)) {
    if (!is.null(cfg$allow_null) && cfg$allow_null) {
      warning(sprintf("Tag not found in DP for %s, returning NULL", name))
      return(NULL)
    } else {
      stop(sprintf("No tag recognised for %s", name))
    }
  }

  conversion_function <- function(x) x
  if (!is.null(cfg$type)) {
    if (cfg$type == "real") {
      conversion_function <- as.numeric
    } else if (cfg$type == "int") {
      conversion_function <- as.integer
    }
  }

  shaped_data <- NULL
  if (is.null(matching_tag_cfg$dims)) {
    shaped_data <- conversion_function(data)
  } else {
    dims <- dim_vars[unlist(matching_tag_cfg$dims)]
    dim_lengths <- unlist(lapply(dims, length))
    dim_names <- lapply(dims, unlist)
    shaped_data <- array(
      conversion_function(unlist(data)),
      dim = dim_lengths,
      dimnames = unname(dim_names)
    )
  }

  list(data = shaped_data, tag = matching_tag_cfg$tag)
}

get_data_from_tag_cfg <- function(tag_cfg, dim_vars, dp) {
  tag_idx <- which(dp[, 1] == sprintf("<%s>", tag_cfg$tag))
  if (length(tag_idx) == 0) return(NULL)

  # Not a fixed offset between different tags and version
  value_idx <- tag_idx
  while (dp$Description[value_idx] != "<Value>") value_idx <- value_idx + 1

  end_idx <- value_idx
  while (dp$Tag[end_idx] != "<End>") end_idx <- end_idx + 1

  start_row <- value_idx
  # we do not count the row with the <End> tag
  end_row <- end_idx - 1

  # PJNZ files have headers:
  # Tag, Description, Notes, Data, X, X.1, X.2, ...
  # we always want to start at Data
  start_column <- which(colnames(dp) == "Data")

  # some datasets have labels in first row such as years so we allow
  # metadata to specify custom start offset
  if (!is.null(tag_cfg$start_offset)) {
    offset <- tag_cfg$start_offset
    start_offset_row <- if (is.null(offset$row)) 0 else offset$row
    start_offset_column <- if (is.null(offset$column)) 0 else offset$column

    start_row <- start_row + start_offset_row
    start_column <- start_column + start_offset_column
  }

  # if no dims then we have a scalar so column end is 1
  # otherwise if user has specified the column_dims use the product of
  # those dim vars or assume the last dim is the length of the columns
  # needed to parse from the PJNZ
  n_cols <- if (is.null(tag_cfg$dims)) {
    1
  } else {
    if (!is.null(tag_cfg$column_dims)) {
      prod(unlist(lapply(tag_cfg$column_dims, function(var) length(dim_vars[[var]]))))
    } else {
      last_dim <- tag_cfg$dims[[length(tag_cfg$dims)]]
      length(dim_vars[[last_dim]])
    }
  }
  end_column <- start_column + n_cols - 1

  skip_rows <- NULL
  skip_columns <- NULL
  if (!is.null(tag_cfg$skip)) {
    if (!is.null(tag_cfg$skip$rows)) {
      # we do not add to the end row index because the end row is
      # calculated by wherever the <End> tag is which will include
      # the skipped rows
      skip_rows <- unlist(tag_cfg$skip$rows)
    }
    if (!is.null(tag_cfg$skip$columns)) {
      end_column <- end_column + length(tag_cfg$skip$columns)
      skip_columns <- unlist(tag_cfg$skip$columns)
    }
  }
  data_rows <- start_row:end_row
  data_columns <- start_column:end_column
  data_rectangle <- dp[data_rows, data_columns]

  if (!is.null(skip_rows)) {
    data_rectangle <- data_rectangle[-skip_rows, ]
  }
  if (!is.null(skip_columns)) {
    data_rectangle <- data_rectangle[, -skip_columns]
  }

  data_rectangle
}


parse_dp <- function(dp) {
  dim_vars <- get_dim_vars(dp)
  metadata <- c(get_years_cfg(), get_pars_metadata(dim_vars, dp))

  ret <- lapply(names(metadata), function(name) {
    get_data_from_cfg(name, metadata[[name]], dim_vars, dp)
  })
  names(ret) <- names(metadata)
  list(data = ret, dim_vars = dim_vars)
}
