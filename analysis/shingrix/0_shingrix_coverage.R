# load packages ----
library("tidyverse")
library("arrow")
library("here")
library("glue")
library("survival")

# Specify analysis parameters ----

# analysis specific parameters
vaccine_name <- "shingrix"
analysis <- "main" 

# define threshold dates
dob_threshold_date <- as.Date("1958-09-01")
threshold_date <- as.Date("2023-09-01")
index_date <- threshold_date
end_date <- threshold_date + years(3) - days(1)

# run functions ----
source(here::here("analysis", "common_code", "utility.R"))
source(here::here("analysis", "common_code", "descriptives.R"))

# create output directory ----
output_dir <- here("output", "shingrix", "coverage")
fs::dir_create(output_dir)

# read and prepare data ----
df_analysis <- read_feather(here(
  "output",
  "shingrix",
  glue("dataset_shingrix_main.arrow")
    )) |>
  mutate(
    age = as.numeric(age),
    age_fct = as.factor(age),

    # define "eligible" as anyone who would turn 65 by end of follow-up
    eligibility = between(
        date_of_birth,
        as.Date("1958-09-01"),
        as.Date("1961-09-01")
      ),
    month_of_birth = floor_date(date_of_birth, "month"), # this should be an identity because currently in OS birth dates are rounded down to 1st of the month
    month_of_birth_fct = fct_reorder(
      scales::label_date(format = "%b '%y")(month_of_birth),
      month_of_birth
    ), # convert to factor
    month_diff_threshold = as.integer(time_length(
      interval(dob_threshold_date, date_of_birth),
      unit = "months"
    )),
    day_diff_threshold = as.integer(date_of_birth - dob_threshold_date),
    registration_duration = as.integer(threshold_date - reg_start_date + 1L),
    censor_date = pmin(reg_end_date, date_of_death, na.rm = TRUE),
    fup_availability = as.numeric(
      pmin(censor_date, end_date, na.rm = TRUE) - threshold_date + 1L
    )
  )

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Cumulative incidence of vaccination since index date ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Cumulative uptake of Shingrix dose 1 by age in months up to one year after the threshold date
cumulative_events(
  df_analysis |> filter(between(age, 64, 65)),
  "month_of_birth_fct",
  precision = 7,
  start_date = threshold_date,
  index_date = threshold_date,
  as.integer(365.25 * 3),
  "shingrix_date_1",
  "Shingrix dose 1 by month"
)

# Cumulative uptake of Shingrix dose 1 by eligibility (born either side of threshold date)
cumulative_events(
  df_analysis,
  "eligibility",
  precision = 7,
  start_date = threshold_date,
  index_date = threshold_date,
  as.integer(365.25 * 3),
  "shingrix_date_1",
  "Shingrix dose 1"
)

# Cumulative uptake of Shingrix dose 2 by age in months up to one year after the threshold date
cumulative_events(
  df_analysis |> filter(between(age, 64, 65)),
  "month_of_birth_fct",
  precision = 7,
  start_date = threshold_date,
  index_date = threshold_date,
  as.integer(365.25 * 3),
  "shingrix_date_2",
  "Shingrix dose 2 by month"
)

# Cumulative uptake of Shingrix dose 2 by eligibility (born either side of threshold date)
cumulative_events(
  df_analysis,
  "eligibility",
  precision = 7,
  start_date = threshold_date,
  index_date = threshold_date,
  as.integer(365.25 * 3),
  "shingrix_date_2",
  "Shingrix dose 2"
)
