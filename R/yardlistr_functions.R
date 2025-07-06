# universal functions -----------------------------------------------------

#' Return list of unique values without changing their type
#'
#' @param x A vector
#' @returns A sorted vector of unique values
#' @keywords internal
unique_values <- function(x) {
  x |>
    dplyr::distinct() |>
    tibble::deframe() |>
    sort()
}

#' Extract time from lubridate datetime
#'
#' @param time A vector of date-time values
#' @returns A vector of times stored as number of seconds since 00:00:00. See `hms::hms()`
#' @keywords internal
get_time <- function(time = lubridate::now()) {
  time |>
    stringr::str_split(" ") |>
    purrr::map_chr(2) |>
    lubridate::hms()
}

#' Calculate tetrad from a datetime
#'
#' A tetrad is one of four parts of a month. Since most months are not
#' a multiple of a week long (except for February in non-leap years),
#' tetrads do not have the same length in days.
#'
#' @param x Vector of datetimes.
#' @param method Method to calculate tetrad.
#'   * `ebird` (the default): The "ebird" method takes day 1-7 as tetrad 1, 8-14 as tetrad 2, 15-21 as tetrad 3 and >21 as tetrad 4. This method is employed by the Cornell Lab of Ornithology in their bar charts on the eBird platform
#'   * `fourths`: the "fourths" method divides a month into four tetrads that are as close to equal in length as possible. A month with 30 days will consists of tetrads of length 7, 8, 7 and 8 in this order. A month with 28 days will consist of tetrads of length 7
#' @returns Vector of tetrad assignments
#' @export
tetrad <- function(x, method = "ebird") {
  # calculate tetrad according to chosen method
  if (method == "fourths") {
    tetrads <- ceiling(lubridate::day(x) / (lubridate::days_in_month(x) / 4))
  } else {
    tetrads <- ceiling(lubridate::day(x) / 7) |> (\(x) replace(x, x == 5, 4))()
  }
}


# specific functions ------------------------------------------------------

#' Make best guess which file from input folder to use
#'
#' @param dir_dat Path to an existing directory containing eBird data
#' @returns The last csv file in ascending alphanumeric order
#' @keywords internal
select_file <- function(dir_dat) {
  # get all files with csv suffix and return last one
  files <- fs::dir_ls(dir_dat, type = "file", glob = "*.csv")
  return(files[length(files)])
}

#' Extract eBird data for specified location
#'
#' @param file_data Path to csv file containing eBird data
#' @param location A string containing the exact name of the location
#' to extract from the data set
#' @param call The environment from which the function is called
#' @returns Tibble in tidy data format
#' @export
clean_ebird_data <- function(file_data, location, call = rlang::caller_env()) {
  # import latest ebird data
  raw_dat <- file_data |>
    readr::read_csv(show_col_types = FALSE) |>
    suppressWarnings()

  # filter data by location
  dat_location <- raw_dat |> filter(.data$Location == location)

  if (length(dat_location$Location) == 0) {
    cli::cli_abort("Location {.val {location}} not found in eBird data")
  }

  # filter data by location and sort by datetime
  dat_cleaned <- dat_location |>
    dplyr::rename(
      "common" = "Common Name",
      "scientific" = "Scientific Name",
      "taxon" = "Taxonomic Order",
      "complete" = "All Obs Reported"
    ) |>
    filter(
      !stringr::str_detect(.data$scientific, "sp.") & # filter out spuhs
        !stringr::str_detect(.data$scientific, "/") # filter out slashes
    ) |>
    tidyr::drop_na(tidyselect::any_of(c("Date", "Time"))) |>
    tidyr::unite("datetime", "Date":"Time", sep = " ") |>
    dplyr::arrange(.data$datetime) |>
    select(c("common", "scientific", "taxon", "datetime", "complete")) |>
    mutate(
      datetime = .data$datetime |> lubridate::ymd_hms(),
      week = lubridate::isoweek(.data$datetime),
      tetrad = tetrad(.data$datetime),
      month = lubridate::month(.data$datetime),
      scientific = stringr::str_replace(.data$scientific, "^", "("),
      scientific = stringr::str_replace(.data$scientific, "$", ")")
    ) |>
    tidyr::unite("species", c("common", "scientific"), sep = " ") |>
    tidyr::unite("month_tetrad", c("month", "tetrad"), remove = FALSE)
}

#' Check if file exists
#'
#' @param dir Path to directory
#' @returns NULL
#'
#' @keywords internal
check_dir_exists <- function(dir) {
  if (!(fs::dir_exists(dir))) {
    cli::cli_abort(c(
      "x" = "Input directory does not exist: {.file {file}}"
    ))
  } else {
    return(invisible(NULL))
  }
}
