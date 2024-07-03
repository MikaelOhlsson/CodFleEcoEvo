library(tidyverse)


Q <- function(z) pnorm(sqrt(2) * z)


b_cod_flounder <- function(m, pars) {
  v <- pars$sigma^2
  growth <- pars$rho * (pars$theta^2 - (pars$zstar - m)^2 - v)
  fishing <- pars$eta * Q((m - pars$phi) / sqrt(2 * v + pars$tau^2))
  hypoxia <- pars$kappa * Q((m - pars$Z) / sqrt(2 * v + pars$nu^2))
  growth - fishing - hypoxia
}


g_cod_flounder <- function(m, pars) {
  v <- pars$sigma^2
  growth <- 2 * pars$rho * v * (pars$zstar - m) / pars$theta^2
  fishing <- pars$eta*v*exp(-(m-pars$phi)^2/(2*v+pars$tau^2))/sqrt(pi*(2*v+pars$tau^2))
  hypoxia <- pars$kappa*v*exp(-(m-pars$Z)^2/(2*v+pars$nu^2))/sqrt(pi*(2*v+pars$nu^2))
  growth - fishing - hypoxia
}


cutoff <- function(n) {
  ifelse(n < 1, (1 * (n > 0)) * (n * n * n * (10 + n * (-15 + 6 * n))), 1)
}


# This function sets up the dynamical equations at any moment of
# time, given the system state at that moment
# Input:
# - time: Moment of time (this is here for compatibility reasons; the
#         system of equations does not actually depend on time explicitly)
# - state: The dynamical state variables (densities and trait means)
#          packaged into a single vector in the following way:
#          state = c(n_1_1, n_1_2, ..., n_1_L, n_2_1, n_2_2, ... n_2_L,
#                    ..., n_S_L, m_1_1, m_1_2, ..., m_1_L, m_2_1, m_2_2,
#                    ..., m_2_L, ..., m_S_L)
#          where n_i_k is the density of species i in patch k, and m_i_k is
#          the trait mean of species i in patch k
# - pars: List of parameters
# Output:
# - Time derivatives of the state variables coerced into a vector (as
#   in the input variable "state" above)
eqs <- function(time, state, pars){
  S <- pars$S # Number of species
  L <- pars$L # Number of habitat patches
  # Arrange species densities and traits into SxL matrices
  n <- matrix(state[1:(S*L)], S, L)
  m <- matrix(state[(S*L + 1):(2*S*L)], S, L)
  # Reserve memory for arrays
  dndt <- matrix(0, S, L) # Time derivative of abundance i in patch k
  dmdt <- matrix(0, S, L) # Time derivative of mean trait i in patch k
  alpha <- array(0, c(S, S, L))
  beta <- array(0, c(S, S, L))
  alphan <- matrix(0, S, L)
  mign <- matrix(0, S, L)
  summig <- matrix(0, S, L)
  betan <- matrix(0, S, L)
  mnmu <- matrix(0, S, L)
  # Calculate ingredient functions
  for (i in 1:S) {
    for (j in 1:S) {
      alpha[i,j,] <- exp(-(m[i,]-m[j,])^2 / (
        2*pars$sigma[i]^2+2*pars$sigma[j]^2+pars$w^2))*pars$w /
        sqrt(2*pars$sigma[i]^2+2*pars$sigma[j]^2+pars$w^2)
      beta[i,j,] <- exp(-(m[i,]-m[j,])^2 / (
        2*pars$sigma[i]^2+2*pars$sigma[j]^2+pars$w^2)) *
        2*pars$sigma[i]^2*pars$w*(m[j,]-m[i,]) / (
          2*pars$sigma[i]^2+2*pars$sigma[j]^2+pars$w^2)^(3/2)
    }
    mign[i,] <- pars$mig[i,,] %*% n[i,] # sum_l mig_ikl * n_il
    summig[i,] <- rep(1, L) %*% pars$mig[i,,] # sum_l mig_ilk
    mnmu[i,] <- pars$mig[i,,] %*% (n[i,]*m[i,]) # sum_l mig_ikl * m_il * n_il
  }
  for (k in 1:L) {
    alphan[,k] <- alpha[,,k] %*% n[,k] # sum_j alpha_ijk * n_jk
    betan[,k] <- beta[,,k] %*% n[,k] # sum_j beta_ijk * n_jk
  }
  ninv <- 1/n # ninv = 1 / n_ik, such that ninv = 0 when n_ik = 0
  ninv[is.infinite(ninv)] <- 0
  b <- b_cod_flounder(m, pars)
  g <- g_cod_flounder(m, pars)
  # Set up the equations
  # The dndt equations are multiplied by a cutoff function, to make
  # behavior smooth as abundances approach 0. The cutoff function is
  # a smoothed-out step function whose derivative exists everywhere.
  dndt <- (n*(b-alphan) + mign - n*summig) * cutoff(n / 1e-10)
  dmdt <- pars$h2*(g - betan + ninv*(mnmu-m*mign))
  # Return equations by first flattening them back into a single vector
  list(c(as.numeric(dndt), as.numeric(dmdt)))
}


# Solve the dynamical equations with given initial conditions and parameters
# Input:
# - pars: List of parameters
# - ic: Initial conditions with the following structure:
#       c(n_1_1, n_2_1, ..., n_S_1, n_1_2, n_2_2, ..., n_S_2, ..., n_S_L,
#         m_1_1, m_2_1, ..., m_S_1, m_1_2, m_2_2, ..., m_S_2, ..., m_S_L)
#        where n_i_k is the density of species i in patch k, and m_i_k is
#        the trait mean of species i in patch k
# - tseq: Vector of sampling point in time for output - e.g., seq(0, 100, by = 0.1)
# - ...: Additional parameters passed to deSolve::ode
# Output:
# - Tidy tibble with solution
integrate_model <- function(pars, ic, tseq, ...) {
  # Table of column names (for future use):
  cnames <- expand_grid(type = c("n", "m"), patch = 1:pars$L, species = 1:pars$S) |>
    mutate(name = str_c(type, species, patch, sep = "_")) |>
    pull(name) |>
    (\(names) c("time", names))()
  # Solve equations & tidy up the solution:
  deSolve::ode(func = pars$model, y = ic, parms = pars, times = tseq, ...) |>
    as.data.frame() |>
    as_tibble() |>
    rename_with(~cnames) |> # Name columns with naming convention "type_species_patch"
    pivot_longer(cols = !time, names_to = "variable", values_to = "v") |>
    separate_wider_delim(variable, delim = "_", names = c("type", "species", "patch")) |>
    pivot_wider(names_from = type, values_from = v) |>
    mutate(species = as.integer(species)) |>
    mutate(sigma = pars$sigma[species],
           species = case_match(species, 1 ~ "cod", 2 ~ "flounder"))
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
    facet_wrap(~patch, nrow = 1, labeller = label_both) +
    labs(x = "Time") +
    scale_y_continuous(name = "Population density", limits = c(0, NA)) +
    scale_colour_manual(values = c("cornflowerblue", "darkseagreen")) +
    theme_cod()
}


plot_trait <- function(sol) {
  sol |>
    ggplot(aes(x = time)) +
    geom_ribbon(aes(ymin = m - sigma, ymax = m + sigma, fill = species), alpha = 0.15) +
    geom_line(aes(y = m, colour = species)) +
    facet_wrap(~patch, nrow = 1, labeller = label_both) +
    scale_colour_manual(values = c("cornflowerblue", "darkseagreen")) +
    scale_fill_manual(values = c("cornflowerblue", "darkseagreen")) +
    labs(x = "Time", y = "Trait value") +
    theme_cod()
}


plot_all <- function(sol) {
  cowplot::plot_grid(plot_density(sol), plot_trait(sol), ncol = 1)
}



tibble(pars = list(list(
  S     = 2, # Number of species
  L     = 2, # Number of patches
  w     = 1, # Competition width
  h2    = c(0.5, 0), # Heritabilities (S)
  sigma = c(0.5, 0.5), # Species trait standard deviations (S)
  theta = matrix(5, 2, 2), # Widths of intrinsic growth functions (SxL)
  rho   = matrix(5, 2, 2), # Maximum intrinsic growth rates (SxL)
  zstar = matrix(1, 2, 2), # Ideal body sizes (SxL)
  eta   = matrix(c(100, 0, 100, 0), 2, 2), # Fishing intensities (SxL)
  phi   = matrix(1, 2, 2), # Fishing body size thresholds
  tau   = matrix(2, 2, 2), # Fishing intensity transition speeds
  Z     = matrix(1, 2, 2), # Hypoxia body size thresholds
  nu    = matrix(1.5, 2, 2), # Hypoxia intensity transition speeds
  kappa = matrix(1, 2, 2), # Maximum effects of hypoxia
  mig   = 0.8 * array(c(0,0,1,0,1,0,0,0), c(2, 2, 2)), # Migration rates (SxLxL)
  model = eqs))
) |>
  mutate(ninit = list(c(130, 0, 130, 167)), # Initial densities n_1_1, n_2_1, ...
         minit = list(c(3.46, 0, 3.46, 0)), # Initial trait means m_1_1, m_2_1, ...
         ic = map2(ninit, minit, c)) |>
  mutate(tseq = list(seq(0, 20, by = 0.01))) |> # Sampling points in time
  mutate(sol = pmap(list(pars, ic, tseq), integrate_model), method = "bdf",
         .keep = "none") |>
  unnest(sol) |>
  mutate(n = ifelse(n > 0, n, 0)) |>
  plot_all()
