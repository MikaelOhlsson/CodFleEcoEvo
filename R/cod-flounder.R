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
  n <- state[1:2] # Species densities
  m <- state[3:4] # Species trait means
  # Ingredients for effects of competition:
  dm <- outer(m, m, FUN = `-`) # Difference matrix of trait means
  v <- pars$sigma^2 # Trait variances
  sv <- outer(v, v, FUN = `+`) # Sum matrix of trait variances
  w2 <- pars$w^2 # Width of competition kernel
  alpha_otherDim <- pars$alpha0 * # Competition from overlap in other traits
    matrix(c(1, pars$alphaI, pars$alphaI, 1), 2, 2)
  alpha <- alpha_otherDim * exp(-dm^2/(2*(sv+w2))) / sqrt(2*pi*(sv+w2))
  # Competitive effect of cod on flounder is zeroed out for simpler model:
  if (pars$feedback == "no feedback") alpha[2,1] <- 0
  beta <- alpha*v*(-dm) / (sv+w2)
  # Ingredients for population density effects of effective intrinsic growth:
  growth <- pars$rho * (pars$theta^2 - (pars$zstar - m)^2 - v) / pars$theta^2
  fishing <- pars$eta * Q((m - pars$phi) / sqrt(2*v + pars$tau^2))
  hypoxia <- pars$kappa * Q((m - pars$zeta) / sqrt(2*v + pars$nu^2))
  b <- growth - fishing - hypoxia
  # Ingredients for trait effects of effective intrinsic growth:
  growth_evo <- 2 * pars$rho * v * (pars$zstar - m) / pars$theta^2
  fishing_evo <- pars$eta * v * exp(-(m-pars$phi)^2/(2*v + pars$tau^2)) /
    sqrt(pi*(2*v + pars$tau^2))
  hypoxia_evo <- pars$kappa * v * exp(-(m-pars$zeta)^2/(2*v + pars$nu^2)) /
    sqrt(pi*(2*v + pars$nu^2))
  g <- growth_evo - fishing_evo - hypoxia_evo
  # Dynamical equations:
  dndt <- n * (b - drop(alpha %*% n)) # Equations for abundances
  dmdt <- pars$h2 * drop(g - beta %*% n) # Equations for trait means
  list(c(dndt, dmdt)) # Return eqs by first flattening them back into a single vector
}


organize_results <- function(sol, pars) {
  as_tibble(as.data.frame(sol)) |>
    rename_with(~paste0("n_", 1:2), 2:3) |> # Rename density columns: n_1, ..., n_S
    rename_with(~paste0("m_", 1:2), 4:5) |> # Trait mean columns:     m_1, ..., m_S
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


generate_params <- function(w = 1, alpha0 = 1, alphaI = 1, theta = 5, rho = 10,
                            zstar = 1, eta = 0, phi = 1.1, tau = 2, zeta = 1, nu = 1.5,
                            kappa = 1, sigma = c(0.5, 0.5), h2 = c(0.5, 0),
                            feedback = "with feedback", model = eqs) {
  list(
    w      = w, # Competition width
    alpha0 = alpha0, # Baseline competition strength
    alphaI = alphaI, # Reduction of competition due to imperfect overlap
    sigma  = sigma, # Species trait standard deviations
    theta  = theta, # Width of intrinsic growth function
    rho    = rho, # Maximum intrinsic growth rate
    h2     = h2, # Heritability; set 2nd entry to 0 to stop flounder evolution
    zstar  = zstar, # Ideal body size for intrinsic growth
    eta    = eta, # Fishing intensity
    phi    = phi, # Fishing body size threshold
    tau    = tau, # Fishing intensity transition speed
    zeta   = zeta, # Hypoxia body size threshold
    nu     = nu, # Hypoxia intensity transition speed
    kappa  = kappa, # Maximum effect of hypoxia
    feedback = feedback, # Set to "no feedback" for no cod-to-flounder interaction
    model  = model
  )
}



crossing(feedback = c("with feedback", "no feedback"),
         fishing_effort = c(0, 11)) |>
  mutate(pars = map2(feedback, fishing_effort,
                     \(f, e) generate_params(feedback = f, eta = e))) |>
  mutate(ninit = list(c(19.1, 17.2)), # Initial species densities
         minit = list(c(0.97, 0)), # Initial species trait means
         ic = map2(ninit, minit, c),
         tseq = list(seq(0, 100, by = 0.1))) |> # Sampling points in time
  mutate(sol = pmap(list(pars, ic, tseq), integrate_model)) |>
  select(feedback, fishing_effort, sol) |>
  mutate(plt = map(sol, plot_all)) |>
  mutate(plt = walk(plt, print))
