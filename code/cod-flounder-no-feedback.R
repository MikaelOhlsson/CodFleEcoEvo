library(tidyverse)



# Sigmoid function from 0 at -inf to 1 at +inf, with inflection point = 1/2 at 0:
Q <- function(z) pnorm(sqrt(2) * z)


eff_intr_growth <- function(m, pars) {
  v <- pars$sigma^2 # Trait variances
  # Components of effective intrinsic growth:
  growth <- pars$rho * (pars$theta^2 - (pars$zstar - m)^2 - v[1]) / pars$theta^2
  fishing <- pars$eta * Q((m - pars$phi) / sqrt(2 * v[1] + pars$tau^2))
  hypoxia <- pars$kappa * Q((m - pars$zeta) / sqrt(2 * v[1] + pars$nu^2))
  growth - fishing - hypoxia # Return total effective intrinsic growth
}


eff_evo_growth <- function(m, pars) {
  v <- pars$sigma^2 # Trait variances
  # Components of trait effects of (effective) intrinsic growth:
  growth_evo <- 2 * pars$rho * v[1] * (pars$zstar - m) / pars$theta^2
  fishing_evo <- pars$eta * v[1] * exp(-(m-pars$phi)^2/(2*v[1]+pars$tau^2)) /
    sqrt(pi*(2*v[1]+pars$tau^2))
  hypoxia_evo <- pars$kappa * v[1] * exp(-(m-pars$zeta)^2/(2*v[1]+pars$nu^2)) /
    sqrt(pi*(2*v[1]+pars$nu^2))
  growth_evo - fishing_evo - hypoxia_evo # Return total trait effects from growth
}


flounder_eqb_density <- function(pars) {
  v <- pars$sigma^2 # Trait variances
  mF <- 0 # Flounder trait mean
  # Components of effective intrinsic growth:
  growth <- pars$rho * (pars$theta^2 - (pars$zstar - mF)^2 - v[2]) / pars$theta^2
  fishing <- pars$eta * Q((mF - pars$phi) / sqrt(2*v[2] + pars$tau^2))
  hypoxia <- pars$kappa * Q((mF - pars$zeta) / sqrt(2*v[2] + pars$nu^2))
  b <- growth - fishing - hypoxia # Intrinsic growth
  aFF <- pars$alpha0 / sqrt(2*pi*(2*v[2] + pars$w^2)) # Self-regulation
  b / aFF # Return equilibrium density
}


dndt <- function(n, m, pars) {
  v <- pars$sigma^2 # Trait variances
  b <- eff_intr_growth(m, pars) # Total effective intrinsic growth
  # Competition coefficients:
  aCC <- pars$alpha0 / sqrt(2*pi*(2*v[1]+pars$w^2))
  aCF <- pars$alpha0 * pars$alphaI * exp(-m^2/(2*(sum(v)+pars$w^2))) /
    sqrt(2*pi*(sum(v)+pars$w^2))
  # Competitor population density:
  Fl <- flounder_eqb_density(pars)
  # Return population's rate of change:
  n * (b - aCC*n - aCF*Fl)
}


dmdt <- function(m, pars) {
  v <- pars$sigma^2 # Trait variances
  g <- eff_evo_growth(m, pars) # Total trait effects from growth
  # Effect of competition on evolution:
  aCF <- pars$alpha0 * pars$alphaI * exp(-m^2/(2*(sum(v)+pars$w^2))) /
    sqrt(2*pi*(sum(v)+pars$w^2))
  bCF <- aCF * (-m) * v[1] / (sum(v)+pars$w^2)
  # Competitor population density:
  Fl <- flounder_eqb_density(pars)
  # Return trait's rate of change:
  pars$h2 * (g - bCF * Fl)
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


plot_density <- function(sol, col = "cornflowerblue") {
  sol |>
    ggplot(aes(x = time, y = n)) +
    geom_line(colour = col) +
    scale_y_continuous(name = "density", limits = c(0, NA)) +
    theme_cod()
}


plot_trait <- function(sol, col = "cornflowerblue") {
  sol |>
    ggplot(aes(x = time)) +
    geom_ribbon(aes(ymin = m - sigma, ymax = m + sigma), fill = col, alpha = 0.15) +
    geom_line(aes(y = m), colour = col) +
    labs(y = "trait") +
    theme_cod()
}


plot_all <- function(sol) {
  cowplot::plot_grid(plot_density(sol), plot_trait(sol), ncol = 1)
}


plot_phase <- function(traits, pars, col = "cornflowerblue") {
  tibble(m = traits) |>
    mutate(dmdt = dmdt(m, pars)) |>
    ggplot(aes(x = m, y = dmdt)) +
    geom_hline(yintercept = 0, alpha = 0.4, linetype = "dashed") +
    geom_line(colour = col) +
    labs(x = expression(paste("Body size (", mu, ")")),
         y = expression(paste("Body size rate of change (", d~mu/d~t, ")"))) +
    theme_cod()
}



tibble(pars = list(list(
  w      = 1, # Competition width
  alpha0 = 1, # Maximum competition between two phenotypes
  alphaI = 1, # Reduction of flounder-to-cod competition from imperfect overlap
  sigma  = c(0.5, 0.5), # Species trait standard deviations
  theta  = 5, # Width of intrinsic growth function
  rho    = 10, # Maximum intrinsic growth rate
  h2     = 0.5, # Heritability; set 2nd entry to 0 to stop flounder evolution
  zstar  = 1, # Ideal body size
  eta    = 12, # Fishing intensity
  phi    = 1.1, # Fishing body size threshold
  tau    = 2, # Fishing intensity transition width
  zeta   = 1, # Hypoxia body size threshold
  nu     = 1.5, # Hypoxia intensity transition width
  kappa  = 1, # Maximum effect of hypoxia
  model  = eqs))
) |>
  mutate(phase_plot = map(pars, \(x) plot_phase(traits = seq(-4, 4, l = 201), x))) |>
  mutate(ninit = 15, # Initial species densities
         minit = 2.5, # Initial species trait means
         ic = map2(ninit, minit, c)) |>
  mutate(tseq = list(seq(0, 100, by = 0.1))) |> # Sampling points in time
  mutate(sol = pmap(list(pars, ic, tseq), integrate_model)) |>
  mutate(dynamics_plot = map(sol, plot_all)) |>
  mutate(plots = map2(dynamics_plot, phase_plot, \(x, y) cowplot::plot_grid(x, y))) |>
  pull(plots)
