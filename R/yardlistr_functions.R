# universal functions -----------------------------------------------------

# return list of unique values without changing their type
unique_values <- function(x) x |> unique() |> tibble::deframe() |> sort()

# extract time from lubridate datetime
get_time <- function(time = lubridate::now()) {
  time %>%
    stringr::str_split(" ") %>%
    purrr::map_chr(2) %>%
    lubridate::hms()
}

# calculate the tetrad (one of four parts of a month) from a datetime
tetrad <- function(x, method = "ebird") {
  # calculate tetrad according to chosen method
  if (method == "fourths") {
    # the "fourths" method divides a month into four tetrads that are as close
    # to equal in length as possible. A month with 30 days will consists of
    # tetrads of length 7, 8, 7 and 8 in this order. A month with 28 days will
    # consist of tetrads of length 7
    tetrads <- ceiling(lubridate::day(x) / (lubridate::days_in_month(x) / 4))
  } else {
    # the "ebird" method takes day 1-7 as tetrad 1, 8-14 as tetrad 2, 15-21 as
    # tetrad 3 and >21 as tetrad 4. This method is employed by the Cornell Lab
    # of Ornithology in their bar charts on the eBird platform
    tetrads <- ceiling(lubridate::day(x) / 7) %>% replace(. == 5, 4)
  }
}


# specific functions ------------------------------------------------------

# make best guess which file from input folder to use
select_file <- function(dir_dat) {
  # get all files with csv suffix and return last one
  files <- fs::dir_ls(dir_dat, type = "file", glob = "*.csv")
  return(files[length(files)])
}

clean_ebird_data <- function(data_file, yard_location,
                             call = rlang::caller_env()) {
  message_file <- stringr::str_glue("Importing Data from ", data_file)
  rlang::inform(message_file)

  # import latest ebird data
  raw_dat <- data_file |>
    readr::read_csv(show_col_types = FALSE) |>
    suppressWarnings()

  # filter data by location
  dat_location <- check_location(raw_dat, yard_location, call)

  # filter data by location and sort by datetime
  dat_cleaned <- dat_location |>
    rename(
      "common" = "Common Name",
      "scientific" = "Scientific Name",
      "taxon" = "Taxonomic Order",
      "complete" = "All Obs Reported"
    ) |>
    filter(
      !stringr::str_detect(scientific, "sp.") & # filter out spuhs
      !stringr::str_detect(scientific, "/") # filter out slashes
    ) |>
    tidyr::drop_na(any_of(c("Date", "Time"))) |>
    tidyr::unite("datetime", Date:Time, sep = " ") |>
    arrange(datetime) |>
    select(c("common", "scientific", "taxon", "datetime", "complete")) |>
    mutate(
      datetime = datetime |> lubridate::ymd_hms(),
      # week = lubridate::isoweek(datetime),
      tetrad = tetrad(datetime),
      month = lubridate::month(datetime),
      scientific = stringr::str_replace(scientific, "^", "("),
      scientific = stringr::str_replace(scientific, "$", ")")
    ) |>
    tidyr::unite("species", c("common", "scientific"), sep = " ") |>
    tidyr::unite("month_tetrad", c("month", "tetrad"), remove = FALSE)
}

# check if specified location is present in data
check_location <- function(raw_dat, yard_location, call) {
  dat_location <- raw_dat |> dplyr::filter(Location == yard_location)

  if (dim(dat_location)[1] == 0) {
    message_location <- stringr::str_glue(
      "Location {yard_location} not found in eBird data"
    )
    rlang::abort(message_location, call = call)
  }

  return(dat_location)
}
