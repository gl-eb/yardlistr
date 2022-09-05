
# yardlistr

<!-- badges: start -->
<!-- badges: end -->

The yardlistr package reads your eBird data and produces a number of different 
plots for a specified location, e.g. your yard (hence the package name)

**Warning**: This package is under development and is provided as is without any
guarantees.
It lacks flexibility and robustenss in various aspects.
If you have questions, issues or suggestions, feel free to open a GitHub issue
or pull request.


## Installation

You can install the development version of yardlistr from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("gl-eb/yardlistr")
```

## How to Use

### File  Organization

First of all, you will need to export your eBird data from
[here](https://ebird.org/downloadMyData).
You will receive an email with a zip file.
When extracted, your data will be in comma-separated format (.csv).
Yardlistr will always choose the last .csv file in the data directory in
ascending alphanumeric order, allowing you to have multiple files.

This example illustrates file selection.
Say you have the following three files in your data directory:

```
ebird_20220802.csv
ebird_20220601.csv
ebird_20210416.csv
```
In ascending alphanumeric order, file ```ebird_20220802.csv``` is last and will
thus be imported by yardlistr.

### Plotting Data

Once your folder structure is set up, you can run yardlistr.
The location name must be exactly the same as on eBird.

``` r
library(yardlistr)
yardlistr("My Location", "path/to/data_directory", "path/to/output_directory")
```
