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


survreg2 <- function(survobj, dist, shift = FALSE) {
  fit <- survreg(survobj ~ 1, dist = dist)
  mu <- fit$coefficients
  sigma <- fit$scale
  alpha <- 1 / sigma
  lambda <- exp(-mu / sigma)

  x_max <- max(as.numeric(survobj))
  x <- 1:x_max

  if (dist %in% c("exponential", "weibull")) {
    surv <- exp(-lambda * x^alpha)
  } else if (dist == "loglogistic") {
    surv <- (1 + lambda * x^alpha)^-1
  }

  if (shift) {
    x <- c(1:shift, x + shift)
    surv <- c(rep(1, shift), surv)
  }

  list("time" = x, "surv" = surv)
}