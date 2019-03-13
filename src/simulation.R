rllogis <- function(n, mu, sigma) {
  # random generation for the log-logistic distribution
  exp(mu + sigma * rlogis(n))
}


simulation <- function(object, sim_duration, demand_rate, n_tables, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (object$dist == "loglogistic") {
    mu <- object$params$mu
    sigma <- object$params$sigma
    surv_time_gen <- function() {round(rllogis(1, mu, sigma)) + object$shift}
  }
  
  hit <- rep(0, sim_duration)
  miss <- integer(sim_duration)
  
  for (t in seq_len(sim_duration)) {
    n_demand <- rpois(1, demand_rate)
    n_enter <- min(n_demand, n_tables - hit[t])
    miss[t] <- n_demand - n_enter
    
    for (i in seq_len(n_enter)) {
      surv_time <- surv_time_gen()
      d <- t:min(t + surv_time - 1, sim_duration)
      hit[d] <- hit[d] + 1
    }
  }
  
  list("hit" = hit, "miss" = miss)
}
