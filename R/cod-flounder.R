library(tidyverse)



# Sigmoid function from 0 at -inf to 1 at +inf, with inflection point = 1/2 at 0:
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


flounder_eqb_density_no_feedback <- function(pars) {
  v <- pars$sigma^2 # Trait variances
  mF <- 0 # Flounder trait mean
  b <- eff_intr_growth(mF, pars)[2] # Intrinsic growth of flounder only
  aFF <- competition(c(mF, mF), pars)$alpha[2,2] # Self-competition of flounder
  b / aFF # Return equilibrium density
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
  b <- eff_intr_growth(m, pars) # Effective intrinsic growth
  g <- eff_evo_growth(m, pars) # Trait effects of (effective) intrinsic growth
  coeff <- competition(m, pars) # Ingredients for effects of competition
  # Dynamical equations:
  dndt <- n * (b - drop(coeff$alpha %*% n)) # Equations for abundances
  dmdt <- pars$h2 * drop(g - coeff$beta %*% n) # Equations for trait means
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


integrate_scenario <- function(pars, ic, tseq, ...) {
  deSolve::ode(func = pars$model, y = ic, parms = pars, times = tseq, ...) |>
    organize_results(pars)
}


# Solve ODEs and put results in a tidy table:
integrate_model <- function(pars, ic, tseq, ...) {
  pars0 <- pars # Define parameter list that is just like the one passed in ...
  pars0$eta <- 0 # ... but with fishing set to zero
  # Solve equations with zero fishing:
  sol0 <- integrate_scenario(pars0, ic, tseq, ...) |>
    mutate(regime = "before fishing")
  # Scrape new initial conditions, as the final state of the no-fishing solution:
  ninit_new <- sol0 |> filter(time == max(tseq)) |> pull(n)
  minit_new <- sol0 |> filter(time == max(tseq)) |> pull(m)
  ic_new <- c(ninit_new, minit_new)
  # Now solve from those initial conditions with fishing turned on:
  sol <- integrate_scenario(pars, ic_new, tseq) |>
    mutate(regime = "with fishing") |>
    mutate(time = time + max(tseq)) # Times are after no-fishing has ended
  sol0 |>
    filter(time < max(time)) |> # Remove duplicated time point
    bind_rows(sol) # Join two tables
}


compute_phase_no_feedback <- function(pars) {
  nF <- flounder_eqb_density_no_feedback(pars) # Equilibrium flounder density
  pars_no_feedback <- pars # Copy parameter list ...
  pars_no_feedback$feedback <- "no feedback" # ... but make sure there's no feedback
  # Generate state vector with cod density at 0 (its value doesn't matter), flounder
  # density at nF, cod trait mean at the input m, and flounder trait mean at 0:
  state <- function(m) c(0, nF, m, 0)
  # Table of possible mean cod trait values:
  tibble(m = seq(-4, 4, l = 201)) |>
    # Compute rate of trait change at each of these possible mean trait values
    # ([[1]] undoes the list; [3] select the cod trait mean equation):
    mutate(dmdt = map_dbl(m, \(m) eqs(0, state(m), pars_no_feedback)[[1]][3]))
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
    geom_vline(xintercept = max(sol$time) / 2, alpha = 0.6, linetype = "dotted") +
    theme_cod()
}


plot_trait <- function(sol) {
  sol |>
    ggplot(aes(x = time)) +
    geom_ribbon(aes(ymin = m - sigma, ymax = m + sigma, fill = species), alpha = 0.15) +
    geom_line(aes(y = m, colour = species)) +
    geom_vline(xintercept = max(sol$time) / 2, alpha = 0.6, linetype = "dotted") +
    labs(x = "Time", y = "Trait") +
    scale_colour_manual(values = c("cornflowerblue", "darkseagreen")) +
    scale_fill_manual(values = c("cornflowerblue", "darkseagreen")) +
    theme_cod()
}


plot_phase_no_feedback <- function(pars) {
  compute_phase_no_feedback(pars) |>
    ggplot(aes(x = m, y = dmdt)) +
    geom_hline(yintercept = 0, alpha = 0.4, linetype = "dashed") +
    geom_line(colour = "cornflowerblue") +
    labs(x = expression(paste("Body size (", mu, ")")),
         y = expression(paste("Rate of change, ", d~mu/d~t))) +
    theme_cod()
}


plot_all <- function(sol, pars) {
  pars0 <- pars
  pars0$eta <- 0
  cowplot::plot_grid(
    plot_density(sol),
    plot_phase_no_feedback(pars0),
    plot_trait(sol),
    plot_phase_no_feedback(pars),
    ncol = 2
  )
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



tibble(feedback = c("with feedback", "no feedback"),
       fishing_effort = 11) |>
  mutate(pars = map2(feedback, fishing_effort,
                     \(f, e) generate_params(feedback = f, eta = e, alphaI = 1))) |>
  mutate(ninit = list(c(19.1, 17.2)), # Initial species densities
         minit = list(c(0.97, 0)), # Initial species trait means
         ic = map2(ninit, minit, c),
         tseq = list(seq(0, 200, by = 0.1))) |> # Sampling points in time
  mutate(sol = pmap(list(pars, ic, tseq), integrate_model)) |>
  select(feedback, fishing_effort, pars, sol) |>
  mutate(plt = map2(sol, pars, plot_all)) |>
  mutate(plt = walk(plt, print))
