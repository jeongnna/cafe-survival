library(tidyverse)
library(survival)
library(RColorBrewer)
source("src/multiplot.R")
source("src/survival-models.R")
source("src/survival-utils.R")


# Preprocessing -----------------------------------------------------------

minimum <- function(..., na.rm = FALSE) {
  values <- cbind(...)
  apply(values, 1, function(x) min(x, na.rm = na.rm))
}

cafe_raw <- read_csv("data/cafe.csv")

cafe <- 
  cafe_raw %>% 
  mutate(x = (out_time - in_time) %>% as.numeric(units = "mins"),
         c = (end_of_study - in_time) %>% as.numeric(units = "mins"),
         event_time = minimum(x, c, na.rm = TRUE),
         delta = 1 - is.na(x)) %>% 
  mutate(party = party_size > 1) %>% 
  select(event_time, delta, cafe_type, party, gender)


# Non-parametric estimation ----------------------------------

survobj <- with(cafe, Surv(event_time, event = delta))
km_survfit <- kaplan_meier(survobj ~ 1)
na_survfit <- nelson_aalen(survobj ~ 1)

config <- list("confidence" = TRUE, 
               "method" = c("geom_step", "geom_step"))
survplot(km_survfit, na_survfit, config = config)

# log-rank test
survdiff(survobj ~ cafe_type, data = cafe)
survdiff(survobj ~ party, data = cafe)
survdiff(survobj ~ gender, data = cafe)
gender_not_mixed <- cafe$gender != "mixed"
survdiff(survobj ~ gender, data = cafe, subset = gender_not_mixed)


# Cox proportional hazard model -------------------------------------------

fit <- survfit(surv ~ cafe_type, data = cafe)
summ <- summary(fit)

time <- summ$time
num_risk <- summ$n.risk
num_event <- summ$n.event
hzd <- num_event / num_risk

tibble(time, hzd, type = c(rep(1, 50), rep(2, 31))) %>% 
  mutate(type = factor(type)) %>% 
  ggplot(aes(time, hzd, col = type, group = type)) +
  geom_point() +
  geom_smooth()

# is_ph <- function(formula, data) {
#   fit <- survfit(formula)
#   summ <- summary(fit)
# }

cox_fit <- coxph(surv ~ cafe_type + party + gender, data = cafe)
cox_fit


# Regression models -------------------------------------------------------

# config for survplot
config = list("confidence" = FALSE,
              "methods" = c("geom_step", rep("geom_line", 3)))
loglog_config = list("confidence" = FALSE,
                     "methods" = c("geom_step", "geom_line"),
                     "palette" = brewer.pal(4, "Set1")[c(1, 4)])

# in univ
univ <- cafe$cafe_type == "in_univ"
univ_regmod <- fit_regression_models(cafe, subset = univ, shift = 5)
multiplot(plotlist = univ_regmod$diagnosis, cols = 3)
survplot(survfit_list = univ_regmod$survfits, config = config)
survplot(survfit_list = univ_regmod$survfits[c(1, 4)], config = loglog_config)

# near offices
office <- cafe$cafe_type == "near_offices"
office_regmod <- fit_regression_models(cafe, subset = office, shift = 0)
multiplot(plotlist = office_regmod$diagnosis, cols = 3)
survplot(survfit_list = office_regmod$survfits, config = config)
survplot(survfit_list = office_regmod$survfits[c(1, 4)], config = loglog_config)
