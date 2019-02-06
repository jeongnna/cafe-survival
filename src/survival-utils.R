survplot <- function(..., survfit_list = NULL, config = NULL) {
  survfits <- c(list(...), survfit_list)
  num_survfits <- length(survfits)

  if (is.null(config)) {
    config <- list()
  }
  if (is.null(config$confidence)) {
    config$confidence <- TRUE
  }
  if (is.null(config$methods)) {
    config$methods <- rep("geom_step", num_survfits)
  }
  if (is.null(config$palette)) {
    config$palette <- brewer.pal(num_survfits, "Set1")[1:num_survfits]
  }
  if (is.null(config$labs)) {
    config$labs <- c("title" = "",
    	               "x" = "Time(min)",
    	               "y" = "Survival probability")
  }

  p <- ggplot()

  for (i in seq_along(survfits)) {
    data <- as_tibble(survfits[[i]])
    method <- get(config$methods[i])
    color <- config$palette[i]

    p <- p + method(aes(time, surv), data = data, col = color)
    if (config$confidence == TRUE) {
      p <-
        p +
        method(aes(time, upper), data = data, col = color, lty = 2) +
        method(aes(time, lower), data = data, col = color, lty = 2)
    }
  }

  p +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = config$labs["title"],
         x = config$labs["x"],
         y = config$labs["y"]) +
    theme_bw()
}


diagnose_dist <- function(formula, print_plot = TRUE) {
  fit <- survfit(formula)
  summ <- summary(fit)

  time <- summ$time
  num_risk <- summ$n.risk
  num_event <- summ$n.event
  cumhzd <- cumsum(num_event / num_risk)
  surv <- summ$surv
  tbl <- tibble(time, cumhzd, surv)

  is_exp <-
    tbl %>%
    mutate(x = time, y = log(surv)) %>%
    ggplot(aes(x, y)) +
    geom_line() +
    geom_smooth(method = "lm", se = FALSE) +
    labs(title = "Is Exponential?", x = "t", y = "log(S(t))") +
    theme_bw()

  is_weibull <-
    tbl %>%
    mutate(x = log(time), y = log(cumhzd)) %>%
    ggplot(aes(x, y)) +
    geom_line() +
    geom_smooth(method = "lm", se = FALSE) +
    labs(title = "Is Weibull?", x = "log(t)", y = "log(H(t))") +
    theme_bw()

  is_loglogistic <-
    tbl %>%
    mutate(x = log(time), y = log(exp(cumhzd) - 1)) %>%
    ggplot(aes(x, y)) +
    geom_line() +
    geom_smooth(method = "lm", se = FALSE) +
    labs(title = "Is log-logistic?", x = "log(t)", y = "log(exp(H(t)-1))") +
    theme_bw()
  plots <- list("exponential" = is_exp,
                "weibull" = is_weibull,
                "loglogistic" = is_loglogistic)
  if (print_plot) {
    multiplot(plotlist = plots, cols = 3)
  } else {
    plots
  }
}


fit_regression_models <- function(data, subset = NULL, shift = NULL) {
  if (is.null(subset)) {
    subset <- 1:nrow(data)
  }
  data <- data[subset, ]

  # Kaplan-Meier (benchmark)
  survobj <- with(data, Surv(event_time, event = delta))
  km_survfit <- kaplan_meier(survobj ~ 1)

  # shift event time
  if (!is.null(shift)) {
    data <-
      data %>%
      mutate(event_time = event_time - shift) %>%
      filter(event_time > 0)
  }

  # create survival object
  survobj <- with(data, Surv(event_time, event = delta))

  # diagnose dists
  diagnosis <- diagnose_dist(survobj ~ 1, print_plot = FALSE)

  # regression models
  exp_survfit <- survreg2(survobj, dist = "exponential", shift = shift)
  weibull_survfit <- survreg2(survobj, dist = "weibull", shift = shift)
  loglog_survfit <- survreg2(survobj, dist = "loglogistic", shift = shift)
  survfits <- list("baseline" = km_survfit,
                   "exponential" = exp_survfit,
                   "weibull" = weibull_survfit,
                   "loglogistic" = loglog_survfit)

  list("diagnosis" = diagnosis, "survfits" = survfits)
}
