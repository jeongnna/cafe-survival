library(tidyverse)
library(survival)
library(KMsurv)
library(RColorBrewer)


# Preprocessing -----------------------------------------------------------

minimum <- function(..., na.rm = FALSE) {
  values <- cbind(...)
  apply(values, 1, function(x) min(x, na.rm = na.rm))
}

cafe_raw <- read_csv("./data/cafe.csv")
cafe_raw

cafe <- 
  cafe_raw %>% 
  mutate(x = (out_time - in_time) %>% as.numeric(units = "mins"),
         c = (end_of_study - in_time) %>% as.numeric(units = "mins"),
         event_time = minimum(x, c, na.rm = TRUE),
         delta = 1 - is.na(x)) %>% 
  mutate(party = party_size > 1) %>% 
  select(event_time, delta, cafe_type, party, gender)
cafe


# Non-parametric Estimation ----------------------------------

kaplan_meier <- function(formula) {
  fit <- survfit(formula)
  summ <- summary(fit)
  
  time <- summ$time
  surv <- summ$surv
  upper <- summ$upper
  lower <- summ$lower
  
  list("time" = time, "surv" = surv, "upper" = upper, "lower" = lower)
}

nelson_aalen <- function(formula) {
  fit <- survfit(formula)
  summ <- summary(fit)
  
  time <- summ$time
  num_risk <- summ$n.risk
  num_event <- summ$n.event
  cumhzd <- cumsum(num_event / num_risk)
  surv <- minimum(1, exp(-cumhzd))
  
  std <- sqrt(cumsum(num_event / num_risk^2))
  upper <- minimum(1, exp(-cumhzd + std * qnorm(0.975)))
  lower <- minimum(1, exp(-cumhzd - std * qnorm(0.975)))
  
  list("time" = time, "surv" = surv, "upper" = upper, "lower" = lower)
}

surv <- Surv(time = cafe$event_time, event = cafe$delta)
km_fit <- kaplan_meier(surv ~ 1)
na_fit <- nelson_aalen(surv ~ 1)

survplot <- function(..., config = NULL) {
  nargs <- length(args <- list(...))
  
  if (is.null(config)) {
    config <- list()
  }
  if (is.null(config$conf_int)) {
    config$conf_int <- TRUE
  }
  if (is.null(config$methods)) {
    config$methods <- rep("geom_step", nargs)
  }
  if (is.null(config$palette)) {
    config$palette <- brewer.pal(nargs, "Set1")[1:nargs]
  }
  if (is.null(config$xylabs)) {
    config$xylabs <- c("Time(min)", "Survival probability")
  }
  
  p <- ggplot()
  
  for (i in seq_along(args)) {
    data <- as_tibble(args[[i]])
    method <- get(config$methods[i])
    color <- config$palette[i]
    
    p <- p + method(aes(time, surv), data = data, col = color)
    if (config$conf_int == TRUE) {
      p <- 
        p + 
        method(aes(time, upper), data = data, col = color, lty = 2) +
        method(aes(time, lower), data = data, col = color, lty = 2)
    }
  }
  
  p +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = config$xylabs[1], y = config$xylabs[2]) +
    theme_bw()
}

config <- list("conf_int" = TRUE, 
               "methods" = c("geom_step", "geom_step"))
survplot(km_fit, na_fit, config = config)

# Log-rank test
survdiff(surv ~ cafe_type, data = cafe)
survdiff(surv ~ party, data = cafe)
survdiff(surv ~ gender, data = cafe)
gender_mask <- cafe$gender != "mixed"
survdiff(surv ~ gender, data = cafe, subset = gender_mask)


# Semi-parametric Estimation ----------------------------------------------

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


# Parametric estimation ---------------------------------------------------


weibull_fit <- survreg(surv ~ 1, dist = "weibull")
weibull_fit
summary(weibull_fit)
class(weibull_fit)


fit <- survfit(surv ~ 1, data = cafe)
summ <- summary(fit)

time <- summ$time
num_risk <- summ$n.risk
num_event <- summ$n.event
survival <- summ$surv
hzd <- num_event / num_risk
cumhzd <- cumsum(hzd)

# Is Exponential?
plot(time, log(survival), type = "l", col = "Blue",
     main="Is Exponential?", xlab = "t", ylab = "log(S(t))")
abline(coef(lm(log(survival) ~ time)), col = "red")

# Is Weibull?
plot(log(time), log(cumhzd),type="l", col = "Blue",
     main = "Is Weibull?", xlab = "log(t)", ylab = "log(H)") ;
abline(coef(lm(log(cumhzd) ~ log(time))), col = "red")

# Is log-logistic?
plot(log(time), log(exp(cumhzd)-1), type="l", col="Blue",
     main="Is log-logistic?", xlab="log(t)", ylab="log(exp(H)-1)");
abline(coef(lm(log(exp(cumhzd)-1)~log(time))), col="red")
