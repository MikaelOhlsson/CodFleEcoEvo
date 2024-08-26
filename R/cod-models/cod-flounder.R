library(tidyverse)



# Sigmoid function from 0 at -inf to 1 at +inf, with inflection point = 1/2 at 0:
Q <- function(z) pnorm(sqrt(2) * z)


cutoff <- function(n) {
  ifelse(n < 1, (1 * (n > 0)) * (n * n * n * (10 + n * (-15 + 6 * n))), 1)
}


# Right-hand side of dynamical equations
# Input:
# - time: Moment of time (this is here for compatibility reasons; the
#         system of equations does not actually depend on time explicitly)
# - state: The dynamical state variables (densities and trait means)
#          packaged into a single vector
# - pars: List of model parameters
# Output:
# - Time derivatives of the state variables coerced into a single vector
eqs <- function(time, state, pars) {
  n <- state[1:pars$S] # Species densities
  m <- state[(pars$S+1):(2*pars$S)] # Species trait means
  # Ingredients for effects of competition:
  dm <- outer(m, m, FUN = `-`) # Difference matrix of trait means
  v <- pars$sigma^2 # Trait variances
  sv <- 2*outer(v, v, FUN = `+`) + diag(rep(pars$w^2, 2)) # Effective variance matrix
  alpha_otherDim <- pars$alpha0 * matrix(c(1, pars$alphaI, pars$alphaI, 1), 2, 2)
  alpha <- alpha_otherDim * exp(-dm^2/sv) / sqrt(pi*sv)
  beta <- alpha*2*v*(-dm) / sv
  # Ingredients for population density effects of effective intrinsic growth:
  growth <- pars$rho * (pars$theta^2 - (pars$zstar - m)^2 - v)
  fishing <- pars$eta * Q((m - pars$phi) / sqrt(2 * v + pars$tau^2))
  hypoxia <- pars$kappa * Q((m - pars$zeta) / sqrt(2 * v + pars$nu^2))
  b <- growth - fishing - hypoxia
  # Ingredients for trait effects of effective intrinsic growth:
  growth_evo <- 2 * pars$rho * v * (pars$zstar - m) / pars$theta^2
  fishing_evo <- pars$eta * v * exp(-(m-pars$phi)^2/(2*v+pars$tau^2)) /
    sqrt(pi*(2*v+pars$tau^2))
  hypoxia_evo <- pars$kappa * v * exp(-(m-pars$zeta)^2/(2*v+pars$nu^2)) /
    sqrt(pi*(2*v+pars$nu^2))
  g <- growth_evo - fishing_evo - hypoxia_evo
  # Dynamical equations:
  dndt <- n * (b - drop(alpha %*% n)) # Equations for abundances
  dmdt <- pars$h2 * drop(g - beta %*% n) # Equations for trait means
  list(c(dndt, dmdt)) # Return eqs by first flattening them back into a single vector
}


organize_results <- function(sol, pars) {
  S <- pars$S # Number of species
  as_tibble(as.data.frame(sol)) |>
    rename_with(~paste0("n_", 1:S), 1 + 1:S) |> # Rename density columns: n_1, ..., n_S
    rename_with(~paste0("m_", 1:S), 1 + S + 1:S) |> # Trait mean columns: m_1, ..., m_S
    pivot_longer(cols = !time, names_to = "variable", values_to = "v") |>
    separate(col = variable, into = c("type", "species"), sep = "_", convert = TRUE) |>
    pivot_wider(names_from = "type", values_from = "v") |> # Separate cols for n and m
    mutate(sigma = pars$sigma[species],
           species = case_match(species, 1 ~ "cod", 2 ~ "flounder"))
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
    ggplot(aes(x = time, y = n, colour = species)) +
    geom_line() +
    labs(x = "Time") +
    scale_y_continuous(name = "Density", limits = c(0, NA)) +
    scale_colour_manual(values = c("cornflowerblue", "darkseagreen")) +
    theme_cod()
}


plot_trait <- function(sol) {
  sol |>
    ggplot(aes(x = time)) +
    geom_ribbon(aes(ymin = m - sigma, ymax = m + sigma, fill = species), alpha = 0.15) +
    geom_line(aes(y = m, colour = species)) +
    labs(x = "Time", y = "Trait") +
    scale_colour_manual(values = c("cornflowerblue", "darkseagreen")) +
    scale_fill_manual(values = c("cornflowerblue", "darkseagreen")) +
    theme_cod()
}


plot_all <- function(sol) {
  cowplot::plot_grid(plot_density(sol), plot_trait(sol), ncol = 1)
}



tibble(pars = list(list(
  S      = 2, # Number of species
  w      = 1, # Competition width
  alpha0 = 1, # Baseline competition strength
  alphaI = 0.1, # Reduction of competition due to imperfect overlap
  sigma  = c(0.5, 0.5), # Species trait standard deviations
  theta  = c(5, 5), # Width of intrinsic growth function
  rho    = c(5, 5), # Maximum intrinsic growth rate
  h2     = c(0.5, 0), # Heritability; set 2nd entry to 0 to stop flounder evolution
  zstar  = c(1, 1), # Ideal body size
  eta    = c(100, 100), # Fishing intensity
  phi    = c(1.1, 1.1), # Fishing body size threshold
  tau    = c(2, 2), # Fishing intensity transition speed
  zeta   = c(1, 1), # Hypoxia body size threshold
  nu     = c(1.5, 1.5), # Hypoxia intensity transition speed
  kappa  = c(1, 1), # Maximum effect of hypoxia
  model  = eqs))
) |>
  mutate(ninit = list(c(300, 320)), # Initial species densities
         minit = list(c(3, 0)), # Initial species trait means
         ic = map2(ninit, minit, c),
         tseq = list(seq(0, 100, by = 0.1))) |> # Sampling points in time
  mutate(sol = pmap(list(pars, ic, tseq), integrate_model), .keep = "none") |>
  unnest(sol) |>
  mutate(n = ifelse(n > 0, n, 0)) |>
  plot_all()
