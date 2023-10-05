
# yardlistr

<!-- badges: start -->
<!-- badges: end -->

The yardlistr package reads your eBird data and produces a number of different plots for a specified location,
e.g. your yard (hence the package name)

**Warning**: This package is under development and is provided as is without any guarantees.
It lacks flexibility and robustenss in various aspects.
If you have questions, issues or suggestions, feel free to open a GitHub issue or pull request.

## Installation

You can install the development version of yardlistr from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("gl-eb/yardlistr")
```

## How to Use

### File  Organization

First of all, you will need to export your eBird data from [here](https://ebird.org/downloadMyData).
You will receive an email with a zip file.
When extracted, your data will be in comma-separated format (.csv).
Yardlistr will always choose the last .csv file in the data directory in ascending alphanumeric order,
allowing you to have multiple files.

This example illustrates file selection.
Say you have the following three files in your data directory:

```
ebird_20220802.csv
ebird_20220601.csv
ebird_20210416.csv
```
In ascending alphanumeric order, file ```ebird_20220802.csv``` is last and will thus be imported by yardlistr.

### Plotting Data

Once your folder structure is set up, you can run yardlistr.
The location name must be exactly the same as on eBird.

``` r
library(yardlistr)
yardlistr("My Location", "path/to/data_directory", "path/to/output_directory")
```

You can also set up a list with multiple locations and pass them to yardlistr one by one.

``` r
library(yardlistr)

dir_input <- "./dir_data"
dir_output <- "./dir_plots"

list_locations <- c(
  "My yard",
  "My friend's yard",
  "My local hotspot"
)

for (location in list_locations) {
  yardlistr::yardlistr(location, dir_input, dir_output)
}
```

#### Species frequency

The frequency plot lists all observed species at the specified location in descending order with a bar representing the percentage of complete checklists that a particular species occurs in.
Species only observed in incidental checklists have a frequency of 0.0%.

![Frequency plot](man/figures/location_frequency.png)

#### Heatmap

The heatmap plot is closely modeled after eBird's bar chart.
For a given species, it displays the percentage of checklists (complete and incidental) that species occurs in for each tetrad of the year.
The higher the frequency of a species in a certain tetrad is, the darker the colour will be.
Each month consists of four tetrads that comprise the following days of the month: 1st - 7th, 8th - 14th, 15th - 21st, 22nd - last.

![Heatmap](man/figures/location_heatmap.png)

#### Species per list

The species per list plot showcases the trajectory of number of species on your location life list versus number of checklists at the location.
In some ways, it resembles a rarefaction curve.

![Species per list](man/figures/location_lists.png)

#### Species over time

This plot is similar to the species per list plot but instead shows how the number of species on your location life list developed over time.
Depending on how regularly you observe birds at the location, the two plots can look more or less similar.

![Species over time](man/figures/location_time.png)

#### Time of day

The time of day plot shows the distribution of checklists throughout the day,
rounded to the nearest hour (e.g. both 6:33 and 7:29 round to 7:00).

![Time of day](man/figures/location_lists_time.png)
