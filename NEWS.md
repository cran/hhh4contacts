# hhh4contacts 0.13.4 (2024-10-07)

* The source code is now hosted in a public Git repository at
  <https://codeberg.org/EE-hub/hhh4contacts>.


# hhh4contacts 0.13.3 (2023-11-06)

* The package now formally depends on **sp** (S4 class used by `map`
  data), not only implicitly via **surveillance**.


# hhh4contacts 0.13.2 (2023-10-23)

* Minor documentation improvements and updates.


# hhh4contacts 0.13.1 (2020-03-20)

* Maintenance for new S3 method lookup in R 4.0.0.
  `NextMethod("confint")` in `confint.fitC()` wouldn't find a method if
  `confint.fitC()` was called directly from `summary.fitC()`.


# hhh4contacts 0.13.0 (2017-06-29)

This is the first version published on CRAN. It requires
[**surveillance**](https://CRAN.R-project.org/package=surveillance) >= 1.14.0.
The main change is in the documentation, which became more comprehensive,
especially for the included datasets.

Other changes include:

* Changed the default `agegroups` parameter in `noroBE()`
  (and `grouping` parameter in `contactmatrix()`)
  such that it corresponds to the six age groups used by Meyer and Held (2017).

* The `aggregateC()` function is now exported.


# hhh4contacts 0.12.6 (2017-01-23)

This version has been used for the publication

> Held L, Meyer S and Bracher J (2017). "Probabilistic forecasting in
> infectious disease epidemiology: The 13th Armitage lecture."
> *Statistics in Medicine*,
> [doi: 10.1002/sim.7363](https://doi.org/10.1002/sim.7363).

## New features

* An additional year of norovirus surveillance counts has been added,
  such that `noroBE()` can now extract data from the time range
  2011-w01 to 2016-w30.

* `noroBE(by = "none")` extracts the overall univariate time series.

* `update()` method for `"fitC"` objects, which enables modification of the
  `"hhh4"` model with automatic refitting of the power parameter of the contact
  matrix.

* `plotHHH4_season_groups()` gained an option to disable confidence intervals.

## Bug fixes

* A typical `drop=FALSE` bug prevented `aggregateC()` from aggregating
  a contact matrix to only two groups.


# hhh4contacts 0.12.1 (2016-07-28)

This is the original version published as supplementary data of

> Meyer S and Held L (2017). "Incorporating social contact data in
> spatio-temporal models for infectious disease spread."
> *Biostatistics*, **18**(2), pp. 338-351,
> [doi: 10.1093/biostatistics/kxw051](https://doi.org/10.1093/biostatistics/kxw051).
