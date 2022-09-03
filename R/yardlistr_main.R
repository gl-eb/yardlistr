#' Main yardlistr functions that calls all other helper functions
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
#' @import dplyr
#' @import ggplot2
#' @importFrom fs path
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

  # get all ebird files
  files <- fs::dir_ls(dir_dat)


  # basic data wrangling ----------------------------------------------------

  # import and clean ebird data
  dat <- clean_ebird_data(files[length(files)], location)

  # list of unique checklist
  checklists <- dat |>
    select(datetime, complete) |>
    unique()
  complete <- checklists |>
    filter(complete == 1)
  # tibble to store count by timepoint
  yardcount <- tibble(
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
      filter(datetime == tibble::deframe(checklists[i, "datetime"])) |>
      select(c("species", "taxon")) |>
      unique()
    # check if any of them have not been seen before
    new_species <- current_species |>
      anti_join(yardlist, by = c("species", "taxon"))
    # add new species to life list
    yardlist <- yardlist |>
      bind_rows(new_species)

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
    select(datetime, month_tetrad) |>
    group_by(month_tetrad) |>
    unique() |>
    count() |>
    ungroup()
  checklists_per_timeunit <- checklists_per_timeunit |>
    full_join(
      timeunits |> anti_join(checklists_per_timeunit, by = "month_tetrad"),
      by = c("month_tetrad", "n")
    ) |>
    mutate(frequency = n / max(n, na.rm = TRUE)) |>
    arrange(month_tetrad) |>
    select(c(month_tetrad, frequency)) |>
    tidyr::pivot_wider(names_from = month_tetrad, values_from = frequency) |>
    tibble::add_column(
      species = c("Checklists"),
      taxon = c(0)
    ) |>
    relocate(c("species", "taxon"), .before = everything())

  # count species frequency per tetrad
  frequency_per_timeunit <- yardlist |>
    # select(-month_tetrad) |>
    bind_cols(matrix(NA, nrow = dim(yardlist)[1], ncol = n_timeunits)) |>
    arrange(taxon)

  # rename columns to numbers
  columns_kept <- c("species", "taxon")
  names(frequency_per_timeunit) <- c(columns_kept, timeunits$month_tetrad)

  # loop through all timeunits and tally bird sightings
  for (i in seq_along(timeunits$month_tetrad)) {
    # filter out current timeunits's data
    dat_timeunit <- dat |>
      filter(month_tetrad == tibble::deframe(timeunits[i, "month_tetrad"]))

    # count number of checklists during timeunit
    n_checklists <- dat_timeunit |>
      select(datetime) |>
      unique_values() |>
      length()

    # normalize species frequency to checklist number
    dat_timeunit <- dat_timeunit |>
      add_count(species) |>
      mutate(frequency = n / n_checklists) |>
      select(c(species, taxon, frequency)) |>
      unique()

    # add species not seen in current tetrad
    dat_timeunit <- dat_timeunit |>
      full_join(
        yardlist |> anti_join(dat_timeunit, by = c("species", "taxon")),
        by = c("species", "taxon")
      ) |>
      arrange(taxon)

    # convert NAs to 0 if at least one checklist present in current tetrad
    if (n_checklists > 0) {
      dat_timeunit <- dat_timeunit |>
        mutate(frequency = tidyr::replace_na(frequency, 0))
    }

    frequency_per_timeunit[, i + length(columns_kept)] <- dat_timeunit |>
      select(frequency) |>
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
    arrange(taxon) |>
    mutate(
      species = species |> forcats::as_factor() |> forcats::fct_rev(),
      frequency = as.numeric(frequency)
    ) |>
    left_join(timeunits, by = c("month_tetrad")) |>
    select(-n)


  # frequency ranking -------------------------------------------------------

  # get overall frequency from frequency per timeunit
  # frequency_year <- frequency_per_timeunit |>
  #   pivot_wider(
  #     id_cols = c("species", "taxon"),
  #     names_from = "month_tetrad",
  #     values_from = "frequency"
  #   ) |>
  #   rowwise(c("species", "taxon")) |>
  #   summarise(frequency = sum(c_across("1_1":last_col()), na.rm = TRUE)) |>
  #   # select(c("species", "taxon", "frequency")) |>
  #   mutate(
  #     frequency = (frequency * stats_checklists$max) / stats_checklists$total
  #   ) |>
  #   filter(species != "Checklists")

  # get overall frequency
  frequency_year <- dat |>
    filter(complete == 1) |>
    count(species, taxon) |>
    full_join(yardlist, by = c("species", "taxon")) |>
    tidyr::replace_na(list(n = 0)) |>
    arrange(desc(n), taxon) %>%
    mutate(
      species = species |> forcats::as_factor() |> forcats::fct_rev(),
      frequency = n / complete |>
        select(datetime) |>
        tibble::deframe() |>
        length(),
      frequency
    )


  # plotting parameters -----------------------------------------------------

  # custom color for plot elements
  clr <- viridisLite::mako(1, begin = 0.7)
  # clr <- "#65A630FF" # eBird green

  # define common theme used for all plots
  plot_theme <- theme_light() +
    theme(
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 13),
      plot.title = element_text(size = 15),
      plot.subtitle = element_text(size = 13)
    )

  # dark theme for heatmap plot
  # plot_theme_heatmap_dark <- plot_theme +
  #   theme_minimal() +
  #   theme(
  #     plot.background = element_rect(fill = "grey10"),
  #     text = element_text(color = "grey87"),
  #     panel.border = element_rect(color = "grey20", fill = "#00000000"),
  #     panel.grid = element_line(color = "grey20"),
  #     axis.text = element_text(color = "grey87"),
  #     strip.text = element_text(size = 11, color = "grey87")
  #   )

  # light theme for heatmap plot
  plot_theme_heatmap_light <- plot_theme +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white"),
      # text = element_text(color = "grey87"),
      # panel.border = element_rect(color = "grey20", fill = "#00000000"),
      # panel.grid = element_line(color = "grey20"),
      # axis.text = element_text(color = "grey87"),
      strip.text = element_text(size = 11)
    )

  # height of plot upon export scales with number of species
  plot_height <- max(10, round(dim(yardlist)[1] / 10) * 5)

  # short location name for file names
  location_short <-
    stringr::str_split(location, stringr::boundary("word"))[[1]][1]


  # plotting ----------------------------------------------------------------

  plot_count_time <- yardcount |>
    ggplot(aes(x = datetime, y = species_num)) +
    geom_step(size = 1, color = clr) +
    labs(
      title = paste0("Yard list at ", location),
      x = "Date",
      y = "Species"
    ) +
    plot_theme
  ggsave(
    filename = path(
      dir_img,
      paste0(location_short, "_time"),
      ext = "png"
    ),
    plot = plot_count_time,
    device = agg_png,
    width = 20,
    height = 12.5,
    units = "cm",
    dpi = 300
  )


  plot_count_lists <- yardcount |>
    ggplot(aes(x = observation, y = species_num)) +
    geom_step(size = 1, color = clr) +
    labs(
      title = paste0("Yard list at ", location),
      x = "Checklists",
      y = "Species"
    ) +
    plot_theme
  ggsave(
    filename = path(
      dir_img,
      paste0(location_short, "_lists"),
      ext = "png"
    ),
    plot = plot_count_lists,
    device = agg_png,
    width = 20,
    height = 12.5,
    units = "cm",
    dpi = 300
  )

  plot_time_lists <- checklists |>
    mutate(time = datetime |> get_time()) |> #
    ggplot(aes(x = time)) +
    geom_histogram(bins = 24, fill = clr, color = "grey13") +
    scale_x_time(
      limits = c(lubridate::hours(0), lubridate::hours(23)),
      labels = scales::label_time(format = "%H:%M")
    ) +
    plot_theme +
    theme(
      axis.title.y = element_blank()
    )
  plot_time_lists
  ggsave(
    filename = path(
      dir_img,
      paste0(location_short, "_lists_time"),
      ext = "png"
    ),
    plot = plot_time_lists,
    device = agg_png,
    width = 20,
    height = 12.5,
    units = "cm",
    dpi = 300
  )

  # plot heatmap
  plot_heatmap <- frequency_per_timeunit |>
    mutate(
      month_label = lubridate::month(month, label = TRUE) |>
        forcats::as_factor()
    ) |>
    ggplot(aes(x = tetrad, y = species, fill = frequency)) +
    facet_grid(cols = vars(month_label), scales = "free_x", switch = "x") +
    geom_raster() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    # scale_fill_viridis_c(option = "magma", na.value = "#00000000") +
    scale_fill_viridis_c(
      option = "mako",
      na.value = "#00000000",
      direction = -1
    ) +
    labs(
      x = "Tetrad",
      fill = "Frequency",
      title = "Species occurence by tetrad",
      subtitle = location
    ) +
    plot_theme_heatmap_light +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_blank(),
      axis.title = element_blank(),
      panel.spacing = unit(0, "null"),
      strip.placement = "outside",
      legend.position = "right"
    )
  # plot_heatmap
  ggsave(
    filename = path(
      dir_img,
      paste0(location_short, "_heatmap"),
      ext = "png"
    ),
    plot = plot_heatmap,
    device = agg_png,
    width = 25,
    height = plot_height,
    units = "cm",
    dpi = 300
  )

  # expand axis limits beyond 100% to accomodate text, if necessary
  plot_frequency_xlim <- c(0, max(1, max(frequency_year$frequency + 0.2)))

  # frequency plot
  plot_frequency <- frequency_year |>
    ggplot(aes(x = frequency, y = species)) +
    scale_x_continuous(
      position = "top",
      limits = plot_frequency_xlim,
      breaks = seq(0, 1, 0.25),
      labels = scales::percent
    ) +
    geom_col(fill = clr) +
    geom_text(
      aes(label = scales::label_percent(accuracy = 0.01)(frequency)),
      nudge_x = 0.1,
      color = "grey30"
    ) +
    labs(
      title = "Species frequency in complete checklists",
      subtitle = location
    ) +
    plot_theme +
    theme(
      plot.title = element_text(hjust = 1),
      plot.subtitle = element_text(hjust = 1),
      axis.title = element_blank()
    )
  # plot_frequency
  ggsave(
    filename = path(
      dir_img,
      paste0(location_short, "_frequency"),
      ext = "png"
    ),
    plot = plot_frequency,
    device = agg_png,
    width = 20,
    height = plot_height,
    units = "cm",
    dpi = 300
  )
}
