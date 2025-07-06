#' Main yardlistr functions for reading, analysing and plotting eBird data
#'
#' @param location (character) Name of location as it appears in eBird data
#' @param dir_dat (character) Relative or absolute path to directory containing
#' eBird data
#' @param dir_img (character) Relative or absolute path to directory to which
#' plots should be saved
#'
#' @return A list containing the tibbles used to generate plots
#' @export
#'
#' @examples
#' \dontrun{
#'   yardlistr("My Location", "./Data", "./Plots")
#' }
yardlistr <- function(location, dir_dat, dir_img) {
  # file system -------------------------------------------------------------

  # set data directories
  dir_dat <- dir_dat |> fs::as_fs_path()
  dir_img <- dir_img |> fs::as_fs_path()

  # check file paths
  check_dir_exists(dir_dat)
  if (!fs::dir_exists(dir_img)) {
    fs::dir_create(dir_img)
  }

  # select file for data import
  file_data <- select_file(dir_dat)

  # basic data wrangling ----------------------------------------------------

  # import and clean ebird data
  dat <- clean_ebird_data(file_data, location)

  # inform user about data import and analysis
  cli::cli_alert_success("Imported Data from {.file {file_data}}")
  cli::cli_alert_info("Analyzing data for {.strong {location}}")

  # list of unique checklist
  checklists <- dat |>
    select("datetime", "complete") |>
    dplyr::distinct()
  complete <- checklists |>
    filter(.data$complete == 1)
  # tibble to store count by timepoint
  yardcount <- tibble::tibble(
    observation = seq_len(dim(checklists)[1]),
    datetime = lubridate::floor_date(checklists$datetime, unit = "day"),
    species_num = numeric(dim(checklists)[1]),
    complete = checklists$complete
  )

  # yardlist ----------------------------------------------------------------

  # tibble to store observed species
  yardlist <- tibble::tibble(
    species = character(),
    taxon = numeric()
  )

  # loop through checklists and calculate number of new species
  for (i in seq_len(dim(checklists)[1])) {
    # get list of species' at current timepoint
    current_species <- dat |>
      filter(.data$datetime == tibble::deframe(checklists[i, "datetime"])) |>
      select(c("species", "taxon")) |>
      dplyr::distinct()
    # check if any of them have not been seen before
    new_species <- current_species |>
      dplyr::anti_join(yardlist, by = c("species", "taxon"))
    # add new species to life list
    yardlist <- yardlist |>
      dplyr::bind_rows(new_species)

    # add number of new species to previous yard count
    if (i == 1) {
      yardcount[i, "species_num"] <- dim(new_species)[1]
    } else {
      yardcount[i, "species_num"] <- yardcount[i - 1, "species_num"] +
        dim(new_species)[1]
    }
  }

  # heatmap -----------------------------------------------------------------

  # list of time units
  n_timeunits <- 48
  timeunits <- tibble::tibble(
    month = rep(seq_len(12), each = 4),
    tetrad = rep(seq_len(4), 12),
    n = rep(NA, n_timeunits)
  ) |>
    tidyr::unite("month_tetrad", c("month", "tetrad"), remove = FALSE)

  # count checklists for each timeunit
  checklists_per_timeunit <- dat |>
    select("datetime", "month_tetrad") |>
    dplyr::group_by(.data$month_tetrad) |>
    dplyr::distinct() |>
    dplyr::count() |>
    dplyr::ungroup()
  checklists_per_timeunit <- checklists_per_timeunit |>
    dplyr::full_join(
      dplyr::anti_join(
        x = timeunits,
        y = checklists_per_timeunit,
        by = dplyr::join_by("month_tetrad")
      ),
      by = c("month_tetrad", "n")
    ) |>
    mutate(frequency = .data$n / max(.data$n, na.rm = TRUE)) |>
    dplyr::arrange(.data$month_tetrad) |>
    select(c("month_tetrad", "frequency")) |>
    tidyr::pivot_wider(
      names_from = "month_tetrad",
      values_from = "frequency"
    ) |>
    tibble::add_column(
      species = c("Checklists"),
      taxon = c(0)
    ) |>
    dplyr::relocate(
      c("species", "taxon"),
      .before = tidyselect::everything()
    )

  # count species frequency per tetrad
  frequency_per_timeunit <- yardlist |>
    dplyr::bind_cols(
      matrix(NA, nrow = dim(yardlist)[1], ncol = n_timeunits)
    ) |>
    dplyr::arrange(.data$taxon) |>
    suppressMessages()

  # rename columns to numbers
  columns_kept <- c("species", "taxon")
  names(frequency_per_timeunit) <- c(columns_kept, timeunits$month_tetrad)

  # loop through all timeunits and tally bird sightings
  for (i in seq_along(timeunits$month_tetrad)) {
    # filter out current timeunits's data
    dat_timeunit <- dat |>
      filter(
        .data$month_tetrad == tibble::deframe(timeunits[i, "month_tetrad"])
      )

    # count number of checklists during timeunit
    n_checklists <- dat_timeunit |>
      select("datetime") |>
      unique_values() |>
      length()

    # normalize species frequency to checklist number
    dat_timeunit <- dat_timeunit |>
      dplyr::add_count(.data$species) |>
      mutate(frequency = .data$n / n_checklists) |>
      select(c("species", "taxon", "frequency")) |>
      dplyr::distinct()

    # add species not seen in current tetrad
    dat_timeunit <- dat_timeunit |>
      dplyr::full_join(
        yardlist |> dplyr::anti_join(dat_timeunit, by = c("species", "taxon")),
        by = dplyr::join_by("species", "taxon")
      ) |>
      dplyr::arrange(.data$taxon)

    # convert NAs to 0 if at least one checklist present in current tetrad
    if (n_checklists > 0) {
      dat_timeunit <- dat_timeunit |>
        mutate(frequency = tidyr::replace_na(.data$frequency, 0))
    }

    frequency_per_timeunit[, i + length(columns_kept)] <- dat_timeunit |>
      select("frequency") |>
      tibble::deframe()
  }

  # merge two tibbles
  frequency_per_timeunit <- checklists_per_timeunit |>
    rbind(frequency_per_timeunit) |>
    tidyr::pivot_longer(
      cols = (length(columns_kept) + 1):(n_timeunits + length(columns_kept)),
      names_to = "month_tetrad",
      values_to = "frequency"
    ) |>
    dplyr::arrange(.data$taxon) |>
    mutate(
      species = .data$species |> forcats::as_factor() |> forcats::fct_rev(),
      frequency = as.numeric(.data$frequency)
    ) |>
    dplyr::left_join(timeunits, by = c("month_tetrad")) |>
    select(-c("n"))

  # frequency ranking -------------------------------------------------------

  # get overall frequency
  frequency_year <- dat |>
    filter(complete == 1) |>
    dplyr::count(.data$species, .data$taxon) |>
    dplyr::full_join(yardlist, by = dplyr::join_by("species", "taxon")) |>
    tidyr::replace_na(list(n = 0)) |>
    dplyr::arrange(dplyr::desc(.data$n), .data$taxon) |>
    mutate(
      species = .data$species |> forcats::as_factor() |> forcats::fct_rev(),
      frequency = .data$n /
        complete |>
          select("datetime") |>
          tibble::deframe() |>
          length()
    )

  # plotting parameters -----------------------------------------------------

  # custom color for plot elements
  clr <- viridisLite::mako(1, begin = 0.7)

  # define common theme used for all plots if user has not set global theme
  no_theme_set <- all.equal(ggplot2::theme_get(), ggplot2::theme_grey())
  if (is.logical(no_theme_set) & isTRUE(no_theme_set)) {
    ggplot2::theme_light(13) |> ggplot2::theme_set()
  }

  # get text color from theme
  clr_text <- ggplot2::theme_get() |> purrr::pluck("axis.text", "colour")

  # height of plot upon export scales with number of species
  plot_height <- max(10, round(dim(yardlist)[1] / 10) * 5)

  # short location name for file names
  location_short <- stringr::str_split(location, stringr::boundary("word")) |>
    purrr::pluck(1, 1)

  # plotting ----------------------------------------------------------------

  plot_count_time <- yardcount |>
    ggplot2::ggplot(aes(x = .data$datetime, y = .data$species_num)) +
    ggplot2::geom_step(linewidth = 1, color = clr) +
    ggplot2::scale_y_continuous(
      breaks = glebrt::breaks_limits(yardcount$species_num)
    ) +
    ggplot2::labs(
      title = paste0("Yard list at ", location),
      x = "Date",
      y = "Species"
    ) +
    ggplot2::theme(
      panel.grid.minor = element_blank()
    )
  ggplot2::ggsave(
    filename = fs::path(
      dir_img,
      paste0(location_short, "_time"),
      ext = "png"
    ),
    plot = plot_count_time,
    device = ragg::agg_png,
    width = 20,
    height = 12.5,
    units = "cm",
    dpi = 300
  )

  plot_count_lists <- yardcount |>
    ggplot2::ggplot(aes(x = .data$observation, y = .data$species_num)) +
    ggplot2::geom_step(linewidth = 1, color = clr) +
    ggplot2::scale_x_continuous(
      breaks = glebrt::breaks_limits(c(0, yardcount$observation))
    ) +
    ggplot2::scale_y_continuous(
      breaks = glebrt::breaks_limits(yardcount$species_num)
    ) +
    ggplot2::labs(
      title = paste0("Yard list at ", location),
      x = "Checklists",
      y = "Species"
    ) +
    ggplot2::theme(
      panel.grid.minor = element_blank()
    )
  ggplot2::ggsave(
    filename = fs::path(
      dir_img,
      paste0(location_short, "_lists"),
      ext = "png"
    ),
    plot = plot_count_lists,
    device = ragg::agg_png,
    width = 20,
    height = 12.5,
    units = "cm",
    dpi = 300
  )

  plot_time_of_day <- checklists |>
    mutate(time = get_time(.data$datetime)) |>
    ggplot2::ggplot(aes(x = .data$time)) +
    ggplot2::stat_bin(
      geom = "bar",
      boundary = lubridate::hms("00:00:00"),
      closed = "left",
      binwidth = lubridate::hours(1),
      fill = clr,
      color = "grey13"
    ) +
    ggplot2::scale_x_time(
      limits = c(lubridate::hms("00:00:00"), lubridate::hms("24:00:00")),
      breaks = scales::breaks_width("4 hours"),
      labels = scales::label_time(format = "%H:%M"),
      expand = ggplot2::expansion(mult = 0.02)
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::labs(
      x = "Time of Day",
      y = "Checklists"
    ) +
    ggplot2::theme(
      panel.grid.major.x = element_blank()
    )
  ggplot2::ggsave(
    filename = fs::path(
      dir_img,
      paste0(location_short, "_time-of-day"),
      ext = "png"
    ),
    plot = plot_time_of_day,
    device = ragg::agg_png,
    width = 20,
    height = 12.5,
    units = "cm",
    dpi = 300
  )

  # plot heatmap
  plot_heatmap <- frequency_per_timeunit |>
    mutate(
      month_label = lubridate::month(.data$month, label = TRUE) |>
        forcats::as_factor()
    ) |>
    ggplot2::ggplot(
      aes(x = .data$tetrad, y = .data$species, fill = .data$frequency)
    ) +
    ggplot2::facet_grid(
      cols = ggplot2::vars(.data$month_label),
      scales = "free_x",
      switch = "x"
    ) +
    ggplot2::geom_raster() +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::scale_fill_viridis_c(
      option = "mako",
      na.value = "#00000000",
      direction = -1
    ) +
    ggplot2::labs(
      x = "Tetrad",
      fill = "Frequency",
      title = "Species occurence by tetrad",
      subtitle = location
    ) +
    ggplot2::theme(
      panel.border = element_blank(),
      panel.spacing = grid::unit(0, "null"),
      panel.grid.major = element_blank(),
      panel.grid.minor.x = element_blank(),
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text = ggplot2::element_text(
        size = 11,
        color = clr_text
      ),
      legend.position = "right",
      axis.ticks = element_blank(),
      axis.text = ggplot2::element_text(color = clr_text),
      axis.text.x = element_blank(),
      axis.title = element_blank()
    )
  ggplot2::ggsave(
    filename = fs::path(
      dir_img,
      paste0(location_short, "_heatmap"),
      ext = "png"
    ),
    plot = plot_heatmap,
    device = ragg::agg_png,
    width = 25,
    height = plot_height,
    units = "cm",
    dpi = 300
  )

  # expand axis limits beyond 100% to accomodate text, if necessary
  plot_frequency_xlim <- c(0, max(1, max(frequency_year$frequency + 0.2)))

  # frequency plot
  plot_frequency <- frequency_year |>
    ggplot2::ggplot(aes(x = .data$frequency, y = .data$species)) +
    ggplot2::scale_x_continuous(
      position = "top",
      limits = plot_frequency_xlim,
      breaks = seq(0, 1, 0.25),
      labels = scales::percent,
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::geom_col(fill = clr) +
    ggplot2::geom_text(
      aes(label = scales::label_percent(accuracy = 0.01)(.data$frequency)),
      nudge_x = 0.1,
      color = "grey30"
    ) +
    ggplot2::labs(
      title = "Species frequency in complete checklists",
      subtitle = location
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 1),
      plot.subtitle = ggplot2::element_text(hjust = 1),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.ticks.y = element_blank()
    )
  ggplot2::ggsave(
    filename = fs::path(
      dir_img,
      paste0(location_short, "_frequency"),
      ext = "png"
    ),
    plot = plot_frequency,
    device = ragg::agg_png,
    width = 20,
    height = plot_height,
    units = "cm",
    dpi = 300
  )

  cli::cli_alert_success("Saved plots to {.file {dir_img}}")

  # return tibbles used for plotting
  list_plot_data <- list(
    "yardcount" = yardcount,
    "checklists" = checklists,
    "frequency_per_timeunit" = frequency_per_timeunit,
    "frequency_year" = frequency_year
  )
}
