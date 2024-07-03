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
  sv <- outer(pars$sigma^2, pars$sigma^2, FUN = `+`) # Sum matrix of trait variances
  alpha <- exp(-dm^2/(2*sv+pars$w^2))*pars$w/sqrt(2*sv+pars$w^2) # alpha matrix
  beta <- alpha*2*pars$sigma^2*(-dm)/(2*sv+pars$w^2) # beta matrix
  # Ingredients for population density effects of effective intrinsic growth:
  growth <- pars$rho * (pars$theta^2 - (pars$zstar - m)^2 - v)
  fishing <- pars$eta * Q((m - pars$phi) / sqrt(2 * v + pars$tau^2))
  hypoxia <- pars$kappa * Q((m - pars$Z) / sqrt(2 * v + pars$nu^2))
  b <- growth - fishing - hypoxia
  # Ingredients for trait effects of effective intrinsic growth:
  growth_evo <- 2 * pars$rho * v * (pars$zstar - m) / pars$theta^2
  fishing_evo <- pars$eta * v * exp(-(m-pars$phi)^2/(2*v+pars$tau^2)) /
    sqrt(pi*(2*v+pars$tau^2))
  hypoxia_evo <- pars$kappa * v * exp(-(m-pars$Z)^2/(2*v+pars$nu^2)) /
    sqrt(pi*(2*v+pars$nu^2))
  g <- growth_evo - fishing_evo - hypoxia_evo
  # Dynamical equations:
  dndt <- n * (b - drop(alpha %*% n)) # Equations for abundances
  dmdt <- pars$h2 * drop(g - beta %*% n) # Equations for trait means
  list(c(dndt, dmdt)) # Return eqs by first flattening them back into a single vector
}


organize_results <- function(sol, pars) {
  dat <- as_tibble(as.data.frame(sol)) # Convert to tibble
  S <- (ncol(dat) - 1) / 2 # Number of species (two cols per species, n and m, plus time)
  dat |>
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
    theme(panel.grid = element_blank(), strip.background = element_rect(fill = "white"))
}


plot_density <- function(sol) {
  sol |>
    ggplot(aes(x = time, y = n, colour = species)) +
    geom_line() +
    scale_y_continuous(name = "density", limits = c(0, NA)) +
    scale_colour_manual(values = c("cornflowerblue", "darkseagreen")) +
    theme_cod()
}


plot_trait <- function(sol) {
  sol |>
    ggplot(aes(x = time)) +
    geom_ribbon(aes(ymin = m - sigma, ymax = m + sigma, fill = species), alpha = 0.15) +
    geom_line(aes(y = m, colour = species)) +
    labs(y = "trait") +
    scale_colour_manual(values = c("cornflowerblue", "darkseagreen")) +
    scale_fill_manual(values = c("cornflowerblue", "darkseagreen")) +
    theme_cod()
}


plot_all <- function(sol) {
  cowplot::plot_grid(plot_density(sol), plot_trait(sol), ncol = 1)
}



tibble(pars = list(list(
  S     = 2, # Number of species
  w     = 1, # Competition width
  sigma = c(0.5, 0.5), # Species trait standard deviations
  theta = c(5, 5), # Width of intrinsic growth function
  rho   = c(5, 5), # Maximum intrinsic growth rate
  h2    = c(0.5, 0), # Heritability; set 2nd entry to 0 to stop flounder evolution
  zstar = c(1, 1), # Ideal body size
  eta   = c(100, 100), # Fishing intensity
  phi   = c(2, 0), # Fishing body size threshold
  tau   = c(2, 1), # Fishing intensity transition speed
  Z     = c(1, 1), # Hypoxia body size threshold
  nu    = c(1.5, 1.5), # Hypoxia intensity transition speed
  kappa = c(1, 1), # Maximum effect of hypoxia
  model = eqs))
) |>
  mutate(ninit = list(c(100, 50)), # Initial species densities
         minit = list(c(3, 0)), # Initial species trait means
         ic = map2(ninit, minit, c)) |>
  mutate(tseq = list(seq(0, 100, by = 0.1))) |> # Sampling points in time
  mutate(sol = pmap(list(pars, ic, tseq), integrate_model)) |>
  unnest(sol) |>
  mutate(n = ifelse(n > 0, n, 0)) |>
  plot_all()
