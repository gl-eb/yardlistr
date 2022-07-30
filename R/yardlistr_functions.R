# functions ---------------------------------------------------------------

# return list of unique values without changing their type
unique_values <- function(x) x |> unique() |> deframe() |> sort()

# extract time from lubridate datetime
get_time <- function(time = now()) {
  time %>%
    str_split(" ") %>%
    map_chr(2) %>%
    hms()
}

# calculate the tetrad (one of four parts of a month) from a datetime
tetrad <- function(x, method = "ebird") {
  # coerce input to vector to make sure downstream operations don't fail
  # x <- as.vector(x)

  # calculate tetrad according to chosen method
  if (method == "fourths") {
    # the "fourths" method divides a month into four tetrads that are as close
    # to equal in length as possible. A month with 30 days will consists of
    # tetrads of length 7, 8, 7 and 8 in this order. A month with 28 days will
    # consist of tetrads of length 7
    tetrads <- ceiling(day(x) / (days_in_month(x) / 4))
  } else {
    # the "ebird" method takes day 1-7 as tetrad 1, 8-14 as tetrad 2, 15-21 as
    # tetrad 3 and >21 as tetrad 4. This method is employed by the Cornell Lab
    # of Ornithology in their bar charts on the eBird platform
    tetrads <- ceiling(day(x) / 7) %>% replace(. == 5, 4)
  }
}

clean_ebird_data <- function(data_file, yard_location) {
  # import latest ebird data
  raw_dat <- data_file |>
    read_csv(show_col_types = FALSE) |>
    suppressWarnings()

  # filter data by location and sort by datetime
  dat <- raw_dat |>
    rename(
      "common" = "Common Name",
      "scientific" = "Scientific Name",
      "taxon" = "Taxonomic Order",
      "complete" = "All Obs Reported"
    ) |>
    filter(
      Location == yard_location & # only selected location
        !str_detect(scientific, "sp.") & # filter out spuhs
        !str_detect(scientific, "/") # filter out slashes
    ) |>
    drop_na(any_of(c("Date", "Time"))) |>
    unite("datetime", Date:Time, sep = " ") |>
    arrange(datetime) |>
    select(c("common", "scientific", "taxon", "datetime", "complete")) |>
    mutate(
      datetime = datetime |> ymd_hms(),
      # week = isoweek(datetime),
      tetrad = tetrad(datetime),
      month = month(datetime),
      scientific = str_replace(scientific, "^", "("),
      scientific = str_replace(scientific, "$", ")")
    ) |>
    unite("species", c("common", "scientific"), sep = " ") |>
    unite("month_tetrad", c("month", "tetrad"), remove = FALSE)
}
