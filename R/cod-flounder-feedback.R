library(tidyverse)



# Sigmoid function from 0 at -inf to 1 at +inf, with inflection point = 1/2 at 0
Q <- function(z) pnorm(sqrt(2) * z)


# Compute effective intrinsic growth rates from the mean traits m and the parameter list
eff_intr_growth <- function(m, pars) {
  v <- pars$sigma^2 # Trait variances
  # Components of effective intrinsic growth:
  growth <- pars$rho * (pars$theta^2 - (pars$zstar - m)^2 - v) / pars$theta^2
  fishing <- pars$eta * Q((m - pars$phi) / sqrt(2*v + pars$tau^2))
  hypoxia <- pars$kappa * Q((m - pars$zeta) / sqrt(2*v + pars$nu^2))
  growth - fishing - hypoxia
}


# Compute rate of trait change from growth, from the mean traits m and the parameter list
eff_evo_growth <- function(m, pars) {
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
  b <- eff_intr_growth(m, pars) # Effective intrinsic growth
  g <- eff_evo_growth(m, pars) # Trait effects of (effective) intrinsic growth
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


integrate_scenario <- function(pars, initCond, tseq) {
  ic <- c(initCond$n, initCond$m)
  deSolve::ode(func = pars$model, y = ic, parms = pars, times = tseq, method = "bdf") |>
    organize_results(pars)
}


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


# Solve ODEs and put results in a tidy table
integrate_model <- function(pars, ic, tseq) {
  pars0 <- pars # Define parameter list that is just like the one passed in via pars,
  pars0$eta <- 0 # but with fishing set to zero
  # Solve equations with zero fishing:
  sol0 <- integrate_scenario(pars0, ic, tseq) |> mutate(regime = "before fishing")
  # Scrape new initial conditions, as the final state of the no-fishing solution:
  ninit_new <- sol0 |> filter(time == max(tseq)) |> pull(n)
  minit_new <- sol0 |> filter(time == max(tseq)) |> pull(m)
  ic_new <- tibble(n = ninit_new, m = minit_new)
  # Now solve from those initial conditions with fishing turned on:
  sol1 <- integrate_scenario(pars, ic_new, tseq) |> mutate(regime = "during fishing") |>
    mutate(time = time + max(tseq)) # Times are after no-fishing has ended
  # Scrape new initial conditions, as the after-state of the fishing solution:
  ninit_after <- sol1 |> filter(time == 2*max(tseq)) |> pull(n)
  minit_after <- sol1 |> filter(time == 2*max(tseq)) |> pull(m)
  ic_after <- tibble(n = ninit_after, m = minit_after)
  # Solve equations with zero fishing after fishing scenario:
  sol <- integrate_scenario(pars0, ic_after, tseq) |> mutate(regime = "after fishing") |>
    mutate(time = time + 2*max(tseq)) # Times are after no-fishing has ended
  sol0 |>
    filter(time < max(time)) |> # Remove duplicated time point
    bind_rows(sol1, sol) # Join two tables
}


plot_density <- function(sol) {
  sol |>
    ggplot(aes(x = time, y = n, colour = species)) +
    geom_line() +
    labs(x = "Time") +
    scale_y_continuous(name = "Density", limits = c(0, NA)) +
    scale_colour_manual(values = c("cornflowerblue", "darkseagreen")) +
    geom_vline(xintercept = c(max(sol$time) / 3, 2* max(sol$time) / 3),
               alpha = 0.6, linetype = "dotted") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_rect(fill = "white"),
          legend.position = "none")
}


plot_trait <- function(sol) {
  sol |>
    ggplot(aes(x = time)) +
    geom_ribbon(aes(ymin = m - sigma, ymax = m + sigma, fill = species), alpha = 0.15) +
    geom_line(aes(y = m, colour = species)) +
    geom_vline(xintercept = c(max(sol$time) / 3, 2* max(sol$time) / 3),
               alpha = 0.6, linetype = "dotted") +
    labs(x = "Time", y = "Trait") +
    scale_colour_manual(values = c("cornflowerblue", "darkseagreen")) +
    scale_fill_manual(values = c("cornflowerblue", "darkseagreen")) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_rect(fill = "white"),
          legend.position = "none")
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



bifurc_data_with_feedback <-
  bifurcation_sim(
    tibble(species = c("cod", "flounder"), n = c(10, 10), m = c(1, 0)),
    fishing_efforts = c(seq(0, 3, by = 0.01), seq(3, 0, by = -0.01)[-1]),
    generate_params(feedback = "with feedback")
  ) |>
  mutate(feedback = "with feedback")

bifurc_data_no_feedback <-
  bifurcation_sim(
    tibble(species = c("cod", "flounder"), n = c(10, 10), m = c(1, 0)),
    fishing_efforts = c(seq(0, 3, by = 0.01), seq(3, 0, by = -0.01)[-1]),
    generate_params(feedback = "no feedback")
  ) |>
  mutate(feedback = "no feedback")

bind_rows(bifurc_data_no_feedback, bifurc_data_with_feedback) |>
  filter(species == "cod") |>
  mutate(species = ifelse(feedback == "no feedback", "Cod; no feedback on flounder",
                          "Cod; with feedback on flounder")) |>
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
#ggsave("../fig/fig-SI-bifurcation.pdf", width = 6, height = 2.8)

tibble(feedback = c("with feedback", "no feedback"), fishing_effort = 2.5) |>
  mutate(pars = map2(feedback, fishing_effort,
                     \(f, e) generate_params(feedback = f, eta = e))) |>
  mutate(ninit = list(c(19.1, 17.2)), # Initial species densities
         minit = list(c(0.97, 0)), # Initial species trait means
         ic = map2(ninit, minit, \(n, m) tibble(n = n, m = m)),
         tseq = list(seq(0, 200, by = 0.1))) |> # Sampling points in time
  mutate(sol = pmap(list(pars, ic, tseq), integrate_model)) |>
  select(feedback, fishing_effort, sol) |>
  crossing(plot_type = c("density", "trait")) |>
  mutate(plot_type = ifelse(plot_type=="trait", list(plot_trait), list(plot_density))) |>
  mutate(plot = map2(sol, plot_type, \(s, p) p(s))) |>
  slice(c(1, 3, 2, 4)) |>
  pull(plot) |>
  cowplot::plot_grid(plotlist = _, nrow = 2, ncol = 2)
