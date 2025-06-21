library(tidyverse)



# Sigmoid function from 0 at -inf to 1 at +inf, with inflection point = 1/2 at 0
Q <- function(z) pnorm(sqrt(2) * z)


# Compute effective intrinsic growth rates from the mean traits m and the parameter list
effective_intr_growth <- function(m, pars) {
  v <- pars$sigma^2 # Trait variances
  # Components of effective intrinsic growth:
  growth <- pars$rho * (pars$theta^2 - (pars$zstar - m)^2 - v) / pars$theta^2
  fishing <- pars$eta * Q((m - pars$phi) / sqrt(2*v + pars$tau^2))
  hypoxia <- pars$kappa * Q((m - pars$zeta) / sqrt(2*v + pars$nu^2))
  growth - fishing - hypoxia
}


# Compute rate of trait change from growth, from the mean traits m and the parameter list
effective_evo_growth <- function(m, pars) {
  v <- pars$sigma^2 # Trait variances
  # Components of trait effects of (effective) intrinsic growth:
  growth_evo <- 2 * pars$rho * v * (pars$zstar - m) / pars$theta^2
  fishing_evo <- pars$eta * v * exp(-(m-pars$phi)^2/(2*v + pars$tau^2)) /
    sqrt(pi*(2*v + pars$tau^2))
  hypoxia_evo <- pars$kappa * v * exp(-(m-pars$zeta)^2/(2*v + pars$nu^2)) /
    sqrt(pi*(2*v + pars$nu^2))
  growth_evo - fishing_evo - hypoxia_evo
}


# Competition & trait change matrices, from mean traits m and the parameter list
competition <- function(m, pars) {
  dm <- outer(m, m, FUN = `-`) # Difference matrix of trait means
  v <- pars$sigma^2 # Trait variances
  sv <- outer(v, v, FUN = `+`) # Sum matrix of trait variances
  alpha_otherDim <- pars$alpha0 * # Competition from overlap in other traits
    matrix(c(1, pars$alphaI, pars$alphaI, 1), 2, 2)
  alpha <- alpha_otherDim * exp(-dm^2/(2*(sv+pars$w^2))) / sqrt(2*pi*(sv+pars$w^2))
  # Competitive effect of cod on flounder is zeroed out for simpler model:
  if (pars$feedback == "no feedback") alpha[2,1] <- 0
  beta <- alpha*v*(-dm) / (sv+pars$w^2)
  list(alpha = alpha, beta = beta)
}


# Right-hand side of dynamical equations
eqs <- function(time, state, pars) {
  n <- state[1:2] # Species densities
  m <- state[3:4] # Species trait means
  b <- effective_intr_growth(m, pars) # Effective intrinsic growth
  g <- effective_evo_growth(m, pars) # Trait effects of (effective) intrinsic growth
  coeff <- competition(m, pars) # Ingredients for effects of competition
  # Dynamical equations:
  dndt <- n * (b - drop(coeff$alpha %*% n)) # Equations for abundances
  dmdt <- pars$h2 * drop(g - coeff$beta %*% n) # Equations for trait means
  list(c(dndt, dmdt)) # Return eqs by first flattening them back into a single vector
}


# Put results in a tidy table
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


# Integrate equations, where initial conditions are given in a table formatted just like
# the final output (this makes repeated integration easier)
integrate_scenario <- function(pars, initCond, tseq) {
  ic <- c(initCond$n, initCond$m) # Extract initial densities and traits from table
  deSolve::ode(func = pars$model, y = ic, parms = pars, times = tseq, method = "bdf") |>
    organize_results(pars)
}


# Generate data for a sequence of fishing efforts, always using the final state of the
# previous run as initial conditions for the next
bifurcation_sim <- function(initCond, fishing_efforts, pars,
                             tseq = seq(0, 200, by = 0.1)) {
  state <- initCond |> mutate(fishing_effort = fishing_efforts[1], .before = 1)
  lastState <- state
  for (eta in fishing_efforts) {
    parsCurrent <- pars
    parsCurrent$eta <- eta
    lastState <- integrate_scenario(parsCurrent, lastState, tseq) |>
      mutate(fishing_effort = eta) |>
      filter(time == max(time)) |>
      select(fishing_effort, species, n, m)
    state <- bind_rows(state, lastState)
  }
  state |> slice(-c(1, 2))
}


generate_params <- function(w = 1, alpha0 = 1, alphaI = 0.5, theta = 5, rho = 5,
                            zstar = 1, eta = 0, phi = 1.1, tau = 2, zeta = 1, nu = 1.5,
                            kappa = 1, sigma = c(0.5, 0.5), h2 = c(0.5, 0),
                            feedback = "with feedback", model = eqs) {
  list(
    w        = w, # Competition width
    alpha0   = alpha0, # Baseline competition strength
    alphaI   = alphaI, # Reduction of competition due to imperfect overlap
    sigma    = sigma, # Species trait standard deviations
    theta    = theta, # Width of intrinsic growth function
    rho      = rho, # Maximum intrinsic growth rate
    h2       = h2, # Heritability; set 2nd entry to 0 to stop flounder evolution
    zstar    = zstar, # Ideal body size for intrinsic growth
    eta      = eta, # Fishing intensity
    phi      = phi, # Fishing body size threshold
    tau      = tau, # Fishing intensity transition speed
    zeta     = zeta, # Hypoxia body size threshold
    nu       = nu, # Hypoxia intensity transition speed
    kappa    = kappa, # Maximum effect of hypoxia
    feedback = feedback, # Set to "no feedback" for no cod-to-flounder interaction
    model    = model
  )
}



bifurcation_data <-
  crossing(
    initCond = list(tibble(species = c("cod", "flounder"), n = c(10, 10), m = c(1, 0))),
    fishing_efforts = list(c(seq(0, 3, by = 0.01), seq(3, 0, by = -0.01)[-1])),
    feedback = c("no feedback", "with feedback")
  ) |>
  mutate(pars = map(feedback, \(x) generate_params(feedback = x))) |>
  mutate(sol = pmap(list(initCond, fishing_efforts, pars), bifurcation_sim))


bifurcation_data |>
  unnest(sol) |>
  filter(species == "cod") |>
  mutate(species = ifelse(feedback == "no feedback", "Cod, without feedback on flounder",
                          "Cod, with feedback on flounder")) |>
  (\(x) bind_rows(x, tibble(
    fishing_effort = unique(x$fishing_effort),
    species = "Flounder", n = NA, m = 0, feedback = NA))
  )() |>
  mutate(species = as_factor(species)) |>
  ggplot(aes(x = fishing_effort, y = m, colour = species, linetype = species)) +
  geom_path() +
  scale_colour_manual(name = NULL, values = c("goldenrod", "steelblue", "gray60")) +
  scale_linetype_manual(name = NULL, values = c("solid", "solid", "dashed")) +
  scale_y_continuous(limits = c(-2.2, 2.2)) +
  labs(x = expression(paste("Fishing effort, ", eta)),
       y = expression(paste("Mean trait value, ", mu))) +
  theme_bw() +
  theme(panel.grid = element_blank())
#ggsave("../fig/fig-SI-bifurcation.pdf", width = 6.2, height = 2.8)
