library(tidyverse)
library(deSolve)



# Sigmoid function: 0 at -inf, 1 at +inf, with inflection point = 1/2 at 0
Q <- function(z) pnorm(sqrt(2) * z)


# Compute effective intrinsic growth rates from the mean traits (m)
# and the list of parameter values (pars)
effective_intr_growth <- function(m, pars) {
  v <- pars$sigma^2 # Trait variances
  # Components of effective intrinsic growth:
  growth <- pars$rho*(pars$theta^2 - (pars$zstar - m)^2 - v) / pars$theta^2
  fishing <- pars$eta*Q((m - pars$phi) / sqrt(2*v + pars$tau^2))
  hypoxia <- pars$kappa*Q((m - pars$zeta) / sqrt(2*v + pars$nu^2))
  growth - fishing - hypoxia
}


# Compute rate of trait change from growth, from the mean traits (m)
# and the list of parameter values (pars)
effective_evo_growth <- function(m, pars) {
  v <- pars$sigma^2 # Trait variances
  # Components of trait effects of (effective) intrinsic growth:
  growth_evo <- 2*pars$rho*v*(pars$zstar - m) / pars$theta^2
  fishing_evo <- pars$eta*v*exp(-(m-pars$phi)^2/(2*v + pars$tau^2)) /
    sqrt(pi*(2*v + pars$tau^2))
  hypoxia_evo <- pars$kappa*v*exp(-(m-pars$zeta)^2/(2*v + pars$nu^2)) /
    sqrt(pi*(2*v + pars$nu^2))
  growth_evo - fishing_evo - hypoxia_evo
}


# Competition & trait change matrices, from the mean traits (m)
# and the list of parameter values (pars)
competition <- function(m, pars) {
  dm <- outer(m, m, FUN = `-`) # Difference matrix of trait means
  v <- pars$sigma^2 # Trait variances
  sv <- outer(v, v, FUN = `+`) # Sum matrix of trait variances
  alpha_other_dim <- pars$alpha0 * # Competition from overlap in
    matrix(c(1, pars$alphaI, pars$alphaI, 1), 2, 2) # non-body size traits
  alpha <- alpha_other_dim * exp(-dm^2/(2*(sv + pars$w^2))) /
    sqrt(2*pi*(sv + pars$w^2))
  beta <- alpha*v*(-dm) / (sv + pars$w^2)
  list(alpha = alpha, beta = beta)
}


# Right hand side of dynamical equations
eqs <- function(time, state, pars) {
  n <- state[1:2] # Species densities
  # If no cod-to-flounder feedback, set flounder density to a constant:
  if (pars$feedback == "no feedback") n[2] <- pars$Fpop
  m <- state[3:4] # Species trait means
  b <- effective_intr_growth(m, pars) # (Effective) intrinsic growth
  g <- effective_evo_growth(m, pars) # Trait effects of intrinsic growth
  coeff <- competition(m, pars) # Ingredients for effects of competition
  # Dynamical equations:
  dndt <- n * (b - drop(coeff$alpha %*% n)) # Equations for abundances
  # If no cod-to-flounder feedback, flounder density does not change:
  if (pars$feedback == "no feedback") dndt[2] <- 0
  dmdt <- pars$h2 * drop(g - coeff$beta %*% n) # Equations for trait means
  list(c(dndt, dmdt))
}


# Put integration results in a tidy table
organize_results <- function(sol, pars) {
  as_tibble(as.data.frame(sol)) |>
    # Rename density columns to n_1, n_2; trait mean columns to m_1, m_2:
    rename(n_1 = `1`, n_2 = `2`, m_1 = `3`, m_2 = `4`) |>
    pivot_longer(cols = !time,
                 names_to = "variable",
                 values_to = "value") |>
    separate(col = variable,
             into = c("type", "species"),
             sep = "_",
             convert = TRUE) |>
    # Separate columns for n (densities) and m (trait means):
    pivot_wider(names_from = type, values_from = value) |>
    # Recode species 1 to "cod", and species 2 to "flounder":
    mutate(species = case_match(species, 1 ~ "cod", 2 ~ "flounder"))
}


# Integrate equations, where initial conditions are given in a table
# formatted just like the final output (makes repeated integration easier)
integrate_scenario <- function(pars, init_cond, tseq, model = eqs) {
  # Extract initial densities and traits from initial condition table:
  ic <- c(init_cond$n, init_cond$m)
  # Integrate equations and put results in tidy format:
  ode(y = ic, times = tseq, func = model, parms = pars, method = "bdf") |>
    organize_results(pars)
}


# Generate data for a sequence of fishing efforts, always using
# the final state of the previous run as initial conditions for the next
bifurcation_sim <- function(init_cond, fishing_efforts, pars,
                            tseq = seq(0, 1000, by = 1)) {
  state <- init_cond |> mutate(fishing_effort = fishing_efforts[[1]])
  last_state <- state
  for (eta in fishing_efforts) {
    pars_current <- pars
    pars_current$eta <- eta
    last_state <- integrate_scenario(pars_current, last_state, tseq) |>
      mutate(fishing_effort = ifelse(species == "cod", eta[1], eta[2])) |>
      filter(time == max(time)) |>
      select(fishing_effort, species, n, m)
    state <- bind_rows(state, last_state)
  }
  state |> slice(-c(1, 2))
}


generate_params <- function(h2 = c(0.5, 0), rho = 5, sigma = c(0.5, 0.5),
                            theta = 5, zstar = 1, eta = c(0, 0), phi = 1.1,
                            tau = 2, kappa = 1, zeta = 1, nu = 1.5,
                            Fpop = 10, w = 1, alpha0 = 1, alphaI = 0.5,
                            feedback = "with feedback") {
  list(
    h2       = h2, # Heritabilities; set 2nd entry to 0 for no flounder evo
    rho      = rho,
    sigma    = sigma, # Vector of trait standard deviations for cod & flounder
    theta    = theta,
    zstar    = zstar,
    eta      = eta, # Vector of fishing intensities for cod and flounder
    phi      = phi,
    tau      = tau,
    kappa    = kappa,
    zeta     = zeta,
    nu       = nu,
    Fpop     = Fpop, # Flounder population density in no-feedback model
    w        = w,
    alpha0   = alpha0,
    alphaI   = alphaI,
    feedback = feedback # Set to "no feedback" for no cod-to-flounder effects
  )
}



bifurcation_data <-
  crossing(
    # Initialize very first simulation with these initial conditions:
    init_cond = list(tibble(species = c("cod", "flounder"),
                            n = c(10, 10),
                            m = c(1, 0))),
    # Gradually increase, then decrease back the fishing effort:
    fishing_efforts = list(
      c(seq(0, 10, by = 0.2), seq(10, 0, by = -0.2)[-1]) |>
        (\(x) list(x, x/2))() |>
        list_transpose()
    )) |>
  # Generate model parameters for both model versions:
  mutate(pars = map(fishing_efforts,
                    \(x) generate_params(eta = x, h2 = c(0.5, 0.5)))) |>
  # Repeatedly integrate models while fishing effort is varied:
  mutate(sol = pmap(list(init_cond, fishing_efforts, pars), bifurcation_sim))


bifurcation_data |>
  unnest(sol) |>
  mutate(fishing_effort = ifelse(species == "cod", fishing_effort, 2*fishing_effort)) |>
  mutate(dir = as_factor(c(rep("increasing fishing effort", 101),
                           rep("decreasing fishing effort", 101)))) |>
  ggplot(aes(x = fishing_effort, y = m, colour = species, linetype = dir)) +
  geom_path() +
  scale_colour_manual(
    name = NULL,
    values = c("goldenrod", "steelblue", "gray60")
  ) +
  scale_linetype_manual(name = NULL, values = c("solid", "dashed")) +
  scale_y_continuous(limits = c(-2.2, 2.2)) +
  labs(x = expression(paste("Fishing effort, ", eta)),
       y = expression(paste("Mean trait value, ", mu))) +
  theme_bw() +
  theme(panel.grid = element_blank())
