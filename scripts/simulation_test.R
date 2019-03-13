source("src/simulation.R")
simres <- simulation(office_regmod$survfits$loglogistic, 1000, 0.45, 30, 123)
plot(simres$occ, type = "l")
sum(simres$miss)



process <- function(n_tables) {
  sim_duration <- 6000
  demand_rate <- 20 / 60
  seed <- 123
  simres <- simulation(object, sim_duration, demand_rate, n_tables, seed)
  sum(simres$miss) / sim_duration
}
object <- univ_regmod$survfits$loglogistic
univ_miss_rate <- sapply(10:30, process)
object <- office_regmod$survfits$loglogistic
office_miss_rate <- sapply(10:30, process)

plot(univ_miss_rate * 60)
plot(office_miss_rate * 60)
