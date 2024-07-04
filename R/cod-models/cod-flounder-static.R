library(tidyverse)



# Sigmoid function from 0 at -inf to 1 at +inf, with inflection point = 1/2 at 0:
Q <- function(z) pnorm(sqrt(2) * z)


cutoff <- function(n) {
  ifelse(n < 1, (1 * (n > 0)) * (n * n * n * (10 + n * (-15 + 6 * n))), 1)
}


eff_intr_growth <- function(m, pars) {
  v <- pars$sigma^2 # Trait variances
  # Components of effective intrinsic growth:
  growth <- pars$rho * (pars$theta^2 - (pars$zstar - m)^2 - v[1])
  fishing <- pars$eta * Q((m - pars$phi) / sqrt(2 * v[1] + pars$tau^2))
  hypoxia <- pars$kappa * Q((m - pars$Z) / sqrt(2 * v[1] + pars$nu^2))
  growth - fishing - hypoxia # Return total effective intrinsic growth
}


eff_evo_growth <- function(m, pars) {
  v <- pars$sigma^2 # Trait variances
  # Components of trait effects of (effective) intrinsic growth:
  growth_evo <- 2 * pars$rho * v[1] * (pars$zstar - m) / pars$theta^2
  fishing_evo <- pars$eta * v[1] * exp(-(m-pars$phi)^2/(2*v[1]+pars$tau^2)) /
    sqrt(pi*(2*v[1]+pars$tau^2))
  hypoxia_evo <- pars$kappa * v[1] * exp(-(m-pars$Z)^2/(2*v[1]+pars$nu^2)) /
    sqrt(pi*(2*v[1]+pars$nu^2))
  growth_evo - fishing_evo - hypoxia_evo # Return total trait effects from growth
}


flounder_eqb_density <- function(pars) {
  v <- pars$sigma^2 # Trait variances
  # Components of effective intrinsic growth:
  growth <- pars$rho * (pars$theta^2 - (pars$zstar - 0)^2 - v[1])
  fishing <- pars$eta * Q((0 - pars$phi) / sqrt(2 * v[2] + pars$tau^2))
  hypoxia <- pars$kappa * Q((0 - pars$Z) / sqrt(2 * v[2] + pars$nu^2))
  b <- growth - fishing - hypoxia # Intrinsic growth
  aFF <- pars$alpha0 * pars$w / sqrt(4 * v[2] + pars$w^2) # Self-regulation
  b / aFF # Return equilibrium density
}


dndt <- function(n, m, pars) {
  v <- pars$sigma^2 # Trait variances
  b <- eff_intr_growth(m, pars) # Total effective intrinsic growth
  # Competition coefficients:
  aCC <- pars$alpha0 * pars$w / sqrt(4 * v[1] + pars$w^2)
  aCF <- pars$alpha0 * pars$alphaF * exp(-m^2/(2*sum(v)+pars$w^2)) * pars$w /
    sqrt(2*sum(v)+pars$w^2)
  # Competitor population density:
  F <- flounder_eqb_density(pars)
  # Return population's rate of change:
  n * (b - aCC * n - aCF * F)
}


dmdt <- function(m, pars) {
  v <- pars$sigma^2 # Trait variances
  g <- eff_evo_growth(m, pars) # Total trait effects from growth
  # Effect of competition on evolution:
  bCF <- pars$alpha0 * pars$alphaF * exp(-m^2/(2*sum(v)+pars$w^2)) *
    2*v[1]*pars$w*(-m) / (2*sum(v)+pars$w^2)^(3/2)
  # Competitor population density:
  F <- flounder_eqb_density(pars)
  # Return trait's rate of change:
  pars$h2 * (g - bCF * F)
}


eqs <- function(time, state, pars) {
  n <- state[1] # Species density
  m <- state[2] # Species trait mean
  list(c(dndt(n, m, pars), dmdt(m, pars)))
}


organize_results <- function(sol, pars) {
  as_tibble(as.data.frame(sol)) |>
    rename(n = `1`, m = `2`) |>
    mutate(sigma = pars$sigma[1])
}


# Solve ODEs and put results in a tidy table:
integrate_model <- function(pars, ic, tseq, ...) {
  deSolve::ode(func = pars$model, y = ic, parms = pars, times = tseq, ...) |>
    organize_results(pars)
}


theme_cod <- function(...) {
  theme_bw(...) +
    theme(panel.grid = element_blank(),
          strip.background = element_rect(fill = "white"),
          legend.position = "none")
}


plot_density <- function(sol) {
  sol |>
    ggplot(aes(x = time, y = n)) +
    geom_line(colour = "cornflowerblue") +
    scale_y_continuous(name = "density", limits = c(0, NA)) +
    theme_cod()
}


plot_trait <- function(sol) {
  sol |>
    ggplot(aes(x = time)) +
    geom_ribbon(aes(ymin = m-sigma, ymax = m+sigma), fill="cornflowerblue", alpha=0.15) +
    geom_line(aes(y = m), colour = "cornflowerblue") +
    labs(y = "trait") +
    theme_cod()
}


plot_all <- function(sol) {
  cowplot::plot_grid(plot_density(sol), plot_trait(sol), ncol = 1)
}


plot_phase <- function(traits, pars) {
  tibble(m = traits) |>
    mutate(dmdt = dmdt(m, pars)) |>
    ggplot(aes(x = m, y = dmdt)) +
    geom_hline(yintercept = 0, alpha = 0.4, linetype = "dashed") +
    geom_line(colour = "cornflowerblue") +
    labs(x = expression(paste("Body size (", mu, ")")),
         y = expression(paste("Body size rate of change (", d~mu/d~t, ")"))) +
    theme_cod()
}



tibble(pars = list(list(
  w      = 1, # Competition width
  alpha0 = 1, # Maximum competition between two phenotypes
  alphaF = 0.02, # Reduction of flounder-to-cod competition from imperfect overlap
  sigma  = c(0.5, 0.5), # Species trait standard deviations
  theta  = 5, # Width of intrinsic growth function
  rho    = 5, # Maximum intrinsic growth rate
  h2     = 0.5, # Heritability; set 2nd entry to 0 to stop flounder evolution
  zstar  = 0, # Ideal body size
  eta    = 0, # Fishing intensity
  phi    = 0, # Fishing body size threshold
  tau    = 0.5, # Fishing intensity transition width
  Z      = 1, # Hypoxia body size threshold
  nu     = 1.5, # Hypoxia intensity transition width
  kappa  = 2, # Maximum effect of hypoxia
  model  = eqs))
) |>
  mutate(phase_plot = map(pars, \(x) plot_phase(traits = seq(-4, 4, l = 201), x))) |>
  mutate(ninit = 40, # Initial species densities
         minit = 4.5, # Initial species trait means
         ic = map2(ninit, minit, c)) |>
  mutate(tseq = list(seq(0, 100, by = 0.1))) |> # Sampling points in time
  mutate(sol = pmap(list(pars, ic, tseq), integrate_model)) |>
  mutate(dynamics_plot = map(sol, plot_all)) |>
  mutate(plots = map2(dynamics_plot, phase_plot, \(x, y) cowplot::plot_grid(x, y))) |>
  pull(plots)
