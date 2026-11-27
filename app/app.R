# Interactive explorer for the cod-flounder eco-evolutionary model.
# Run from the repository root with shiny::runApp("app").

library(shiny)
library(bslib)
library(tidyverse)
library(deSolve)


defaults <- list(
  rho = 5, theta = 5, zstar = 1,
  eta = 3, etaF = 0, phi = 1, tau = 0.5,
  kappa = 2, zeta = 1, nu = 1.5,
  alpha0 = 0.707, alphaI = 0.5, w = 0.33,
  sC = 0.5, sF = 0.5, h2 = 0.5,
  NF = 12,
  tmax = 150, fishing = c(0, 50), muC0 = 1.5, NC0 = 5, burnin = TRUE,
  eta_max_bif = 6
)

col_cod <- "cornflowerblue"
col_flounder <- "darkseagreen"
col_fishing <- "firebrick"


# Model ----------------------------------------------------------------------

# Sigmoid: (erf(z) + 1) / 2
Q <- function(z) pnorm(sqrt(2) * z)

# Components of the per capita growth rate (Eq. 3) for a species with mean
# trait m, trait variance v, and fishing effort eta
growth_term <- function(m, v, p) {
  p$rho * (p$theta^2 - (p$zstar - m)^2 - v) / p$theta^2
}
fishing_term <- function(m, v, eta, p) {
  eta * Q((m - p$phi) / sqrt(2*v + p$tau^2))
}
hypoxia_term <- function(m, v, p) {
  p$kappa * Q((m - p$zeta) / sqrt(2*v + p$nu^2))
}

# Corresponding selection terms (inside the brackets of Eq. 4)
growth_sel <- function(m, v, p) 2 * p$rho * v * (p$zstar - m) / p$theta^2
fishing_sel <- function(m, v, eta, p) {
  -eta * v * exp(-(m - p$phi)^2 / (2*v + p$tau^2)) / sqrt(pi*(2*v + p$tau^2))
}
hypoxia_sel <- function(m, v, p) {
  -p$kappa * v * exp(-(m - p$zeta)^2 / (2*v + p$nu^2)) /
    sqrt(pi*(2*v + p$nu^2))
}

# Competition coefficient (Eq. 5); `other` is 1 within species, alphaI between
alpha_fn <- function(mi, mj, vi, vj, other, p) {
  S <- vi + vj + p$w^2
  p$alpha0 * other * exp(-(mi - mj)^2 / (2*S)) / sqrt(2*pi*S)
}

# Selection on species i from competition with species j (density Nj)
comp_sel <- function(mi, mj, vi, vj, Nj, other, p) {
  (mi - mj) * vi / (vi + vj + p$w^2) * alpha_fn(mi, mj, vi, vj, other, p) * Nj
}

# Densities are kept from falling below 1e-100 (effectively extinct),
# by smoothly damping negative growth rates as that floor is approached; this
# has no effect on the dynamics at any density above ~1e-95
log_floor <- log(1e-100)
floor_rate <- function(lN, rate) {
  # (the exponent is capped so that trial steps of the solver far below the
  # floor do not overflow)
  if (rate >= 0) rate else rate * -expm1(min(log_floor - lN, 1))
}

# Per capita growth rates of cod and flounder, and the rate of change of cod
# mean size (vectorized). With dynamic = FALSE, flounder density is p$NF.
rates <- function(lNC, muC, lNF, eta, etaF, p, dynamic) {
  vC <- p$sC^2
  vF <- p$sF^2
  # (exponents are capped so that wild trial steps of the solver cannot
  # overflow to Inf; real densities never come near e^50)
  NC <- exp(pmin(lNC, 50))
  NF <- if (dynamic) exp(pmin(lNF, 50)) else p$NF
  muF <- 0
  aCC <- alpha_fn(muC, muC, vC, vC, 1, p)
  aCF <- alpha_fn(muC, muF, vC, vF, p$alphaI, p)
  bC <- growth_term(muC, vC, p) - fishing_term(muC, vC, eta, p) -
    hypoxia_term(muC, vC, p)
  out <- list(
    rC = bC - aCC * NC - aCF * NF,
    dmuC = p$h2 * (growth_sel(muC, vC, p) + fishing_sel(muC, vC, eta, p) +
                     hypoxia_sel(muC, vC, p) +
                     comp_sel(muC, muF, vC, vF, NF, p$alphaI, p))
  )
  if (!dynamic) return(out)
  aFF <- alpha_fn(muF, muF, vF, vF, 1, p)
  aFC <- alpha_fn(muF, muC, vF, vC, p$alphaI, p)
  bF <- growth_term(muF, vF, p) - fishing_term(muF, vF, etaF, p) -
    hypoxia_term(muF, vF, p)
  c(out, list(rF = bF - aFF * NF - aFC * NC))
}

# Right-hand side of the dynamics. Densities are integrated on a log scale
# (state: lNC = log N_C, muC, and lNF = log N_F if dynamic), which keeps them
# positive and lets the solver cope with populations driven to tiny values.
rhs <- function(t, y, parms) {
  dyn <- parms$dynamic
  r <- rates(y[["lNC"]], y[["muC"]], if (dyn) y[["lNF"]] else NA,
             parms$eta, parms$etaF, parms$p, dyn)
  d <- c(floor_rate(y[["lNC"]], r$rC), r$dmuC)
  if (dyn) d <- c(d, floor_rate(y[["lNF"]], r$rF))
  list(d)
}

# Integrate from y0 = c(NC, muC[, NF]); returns densities on natural scale
integrate <- function(y0, times, eta, etaF, p, dynamic) {
  z0 <- c(lNC = log(y0[["NC"]]), muC = y0[["muC"]])
  if (dynamic) z0 <- c(z0, lNF = log(y0[["NF"]]))
  out <- ode(z0, times, rhs,
             list(p = p, eta = eta, etaF = etaF, dynamic = dynamic),
             method = "lsoda", hmax = 5) |>
    as.data.frame() |>
    as_tibble() |>
    mutate(NC = exp(lNC), .after = time) |>
    select(-lNC)
  if (dynamic) out <- out |> mutate(NF = exp(lNF)) |> select(-lNF)
  out
}

# Simulate with fishing switched on during [t1, t2)
simulate <- function(p, dynamic, t1, t2, tmax, y0, burnin) {
  if (burnin) { # Start from the equilibrium reached without any fishing
    y0 <- unlist(tail(integrate(y0, c(0, 2000), 0, 0, p, dynamic), 1)[-1])
  }
  times <- sort(unique(c(seq(0, tmax, length.out = 601), t1, t2)))
  segments <- list(c(0, t1, 0), c(t1, t2, 1), c(t2, tmax, 0))
  out <- list()
  for (s in segments) {
    tt <- times[times >= s[1] & times <= s[2]]
    if (length(tt) < 2) next
    res <- integrate(y0, tt, s[3] * p$eta, s[3] * p$etaF, p, dynamic)
    y0 <- unlist(tail(res, 1)[-1])
    out <- c(out, list(res))
  }
  res <- bind_rows(out) |> distinct(time, .keep_all = TRUE)
  if (!dynamic) res$NF <- p$NF
  res
}


# Analysis: equilibria, thresholds, bifurcation diagram ----------------------
#
# For a given cod mean size m, the densities settle to their ecological
# equilibrium: with fixed N_F, N_C = (b_C - a_CF N_F) / a_CC; with dynamic
# N_F, (N_C, N_F) solve the 2x2 Lotka-Volterra system. These densities are
# linear in the fishing effort eta, and therefore so is d mu_C / dt evaluated
# at them. Equilibria for a given eta are thus the solutions of eta = H(mu)
# (see code/parameter-sweeps.R), so the bifurcation diagram is the curve
# (H(mu), mu), and the tipping points are the local extrema of H.
#
# In the dynamic model, flounder are fished in proportion to cod (at the ratio
# eta_F / eta set in the sidebar), so that the analysis at the current eta
# matches the fishing period of the simulation.

etaF_at <- function(eta, p) if (p$eta > 0) eta * p$etaF / p$eta else 0 * eta

# Flounder equilibrium density without cod
flounder_alone <- function(eta, p) {
  vF <- p$sF^2
  (growth_term(0, vF, p) - fishing_term(0, vF, etaF_at(eta, p), p) -
     hypoxia_term(0, vF, p)) / alpha_fn(0, 0, vF, vF, 1, p)
}

# Ecological equilibrium densities at cod mean size m (vectorized in m and
# eta). By default, this is the coexistence solution, in which a density can
# be negative (meaning that species cannot persist there). With clamp = TRUE,
# cod that cannot persist are set to zero (with flounder alone at their own
# equilibrium in the dynamic model).
eco_eq <- function(m, eta, p, dynamic, clamp = FALSE) {
  vC <- p$sC^2
  vF <- p$sF^2
  aCC <- alpha_fn(m, m, vC, vC, 1, p)
  aCF <- alpha_fn(m, 0, vC, vF, p$alphaI, p)
  bC <- growth_term(m, vC, p) - fishing_term(m, vC, eta, p) -
    hypoxia_term(m, vC, p)
  if (!dynamic) {
    NF <- p$NF + 0 * bC
    NC <- (bC - aCF * NF) / aCC
    return(list(NC = if (clamp) pmax(NC, 0) else NC, NF = NF))
  }
  aFF <- alpha_fn(0, 0, vF, vF, 1, p)
  aFC <- alpha_fn(0, m, vF, vC, p$alphaI, p)
  bF <- aFF * flounder_alone(eta, p)
  det <- aCC * aFF - aCF * aFC
  det[abs(det) < 1e-12] <- NA # Cod and flounder identical (alphaI = 1, m = 0)
  NC <- (aFF * bC - aCF * bF) / det
  NF <- (aCC * bF - aFC * bC) / det
  if (clamp) {
    absent <- is.na(NC) | NC <= 0
    NF0 <- pmax(bF / aFF, 0) + 0 * NC
    NC[absent] <- 0
    NF[absent] <- NF0[absent]
  }
  list(NC = NC, NF = NF)
}

# d mu_C / dt for given flounder density NF
dmu_NF <- function(m, eta, NF, p) {
  vC <- p$sC^2
  p$h2 * (growth_sel(m, vC, p) + fishing_sel(m, vC, eta, p) +
            hypoxia_sel(m, vC, p) +
            comp_sel(m, 0, vC, p$sF^2, NF, p$alphaI, p))
}

# d mu_C / dt with densities at their ecological equilibrium
dmu_eq <- function(m, eta, p, dynamic, clamp = TRUE) {
  dmu_NF(m, eta, eco_eq(m, eta, p, dynamic, clamp)$NF, p)
}

# Fishing effort at which m is an equilibrium with cod present: eta = H(m)
H_fn <- function(m, p, dynamic) {
  d0 <- dmu_eq(m, 0, p, dynamic, clamp = FALSE)
  d1 <- dmu_eq(m, 1, p, dynamic, clamp = FALSE)
  -d0 / (d1 - d0)
}

# Same, with cod absent (dynamic model; flounder alone at equilibrium). This
# branch is only an equilibrium where cod cannot invade.
H_absent <- function(m, p) {
  d0 <- dmu_NF(m, 0, flounder_alone(0, p), p)
  d1 <- dmu_NF(m, 1, flounder_alone(1, p), p)
  -d0 / (d1 - d0)
}

mu_grid <- function(p) {
  seq(min(-4, p$zstar - 4), max(4, p$zstar + 4), by = 0.002)
}

# Largest real part of the eigenvalues of the Jacobian of the dynamic model,
# at equilibria with both species present (vectorized)
max_eigen <- function(m, NC, NF, eta, p) {
  if (length(m) == 0) return(numeric())
  f <- function(X) {
    r <- rates(X[, 1], X[, 2], X[, 3], eta, etaF_at(eta, p), p, TRUE)
    cbind(r$rC, r$dmuC, r$rF)
  }
  X <- cbind(log(NC), m, log(NF))
  h <- 1e-6
  D <- lapply(1:3, \(i) {
    E <- matrix(0, nrow(X), 3)
    E[, i] <- h
    (f(X + E) - f(X - E)) / (2 * h)
  })
  map_dbl(seq_along(m), \(k) {
    J <- cbind(D[[1]][k, ], D[[2]][k, ], D[[3]][k, ])
    max(Re(eigen(J, only.values = TRUE)$values))
  })
}

# Equilibria of mu_C at fishing effort eta, with their stability and densities
equilibria <- function(eta, p, dynamic) {
  g <- mu_grid(p)
  d <- dmu_eq(g, eta, p, dynamic)
  i <- which(diff(sign(d)) != 0)
  if (length(i) == 0) {
    return(tibble(mu = numeric(), stable = logical(), NC = numeric(),
                  NF = numeric()))
  }
  mu <- map_dbl(i, \(k) uniroot(\(m) dmu_eq(m, eta, p, dynamic),
                                g[c(k, k + 1)], tol = 1e-10)$root)
  N <- eco_eq(mu, eta, p, dynamic, clamp = TRUE)
  stable <- d[i] > 0
  if (dynamic) { # Full Jacobian where both species are present
    both <- N$NC > 0 & N$NF > 0
    stable[both] <- max_eigen(mu[both], N$NC[both], N$NF[both], eta, p) < 0
  }
  tibble(mu = mu, stable = stable, NC = N$NC, NF = N$NF)
}

# The equilibrium curve eta = H(mu) on the grid, prepared for finding the
# tipping points: H is set to NA at poles (where d mu_C / dt does not depend
# on eta, so H jumps between -Inf and Inf) and, in the dynamic model, where
# the coexistence solution has a negative density (not an equilibrium).
# Stability along the curve follows from d(d mu_C / dt) / d mu = -H' B, where
# B = d(d mu_C / dt) / d eta: stable where H' B > 0 (with fixed N_F, B < 0,
# so simply where H' < 0).
eq_curve <- function(p, dynamic) {
  g <- mu_grid(p)
  B <- dmu_eq(g, 1, p, dynamic, clamp = FALSE) -
    dmu_eq(g, 0, p, dynamic, clamp = FALSE)
  h <- H_fn(g, p, dynamic)
  h[!is.finite(h)] <- NA
  pole <- which(diff(sign(B)) != 0)
  h[c(pole, pole + 1)] <- NA
  if (dynamic) {
    N <- eco_eq(g, h, p, dynamic)
    h[is.na(N$NC) | N$NC <= 0 | N$NF <= 0] <- NA
  }
  # Slope of H (one-sided next to gaps in the curve)
  dh <- coalesce((lead(h) - lag(h)) / 2, lead(h) - h, h - lag(h))
  tibble(mu = g, eta = h, stable = !is.na(h) & !is.na(dh) & dh * B > 0)
}

# Tipping points and regime. Starting from the large-cod state without
# fishing, the stable branch is followed towards larger eta: eta_crit is where
# it ends in a fold. The stable state the system then jumps to is followed
# back towards smaller eta: eta_rec is where that branch ends.
fold_points <- function(p, dynamic) {
  d <- eq_curve(p, dynamic)
  g <- d$mu
  h <- d$eta
  eq0 <- equilibria(0, p, dynamic) |> filter(stable)
  mu0 <- if (nrow(eq0) > 0) max(eq0$mu) else NA_real_
  none <- list(eta_crit = NA_real_, eta_rec = NA_real_, mu0 = mu0)
  if (is.na(mu0)) return(c(none, regime = "no stable state without fishing"))
  if (mu0 < 0) return(c(none, regime = "cod small even without fishing"))
  # Stable run of grid points containing a given index
  run <- function(i) {
    a <- i
    while (a > 1 && d$stable[a - 1]) a <- a - 1
    b <- i
    while (b < nrow(d) && d$stable[b + 1]) b <- b + 1
    a:b
  }
  refine <- function(i, maximum) {
    optimize(\(m) H_fn(m, p, dynamic), g[c(i - 1, i + 1)],
             maximum = maximum)$objective
  }
  i0 <- which.min(abs(g - mu0))
  if (!d$stable[i0]) return(c(none, regime = "complex (check bifurcation diagram)"))
  # Large-cod branch: eta increases towards smaller mu, up to its lower end a
  a <- min(run(i0))
  if (a == 1 || is.na(h[a - 1])) {
    if (dynamic && a > 1 && !is.na(H_fn(g[a - 1], p, dynamic))) {
      return(c(none, regime = "large cod die out before tipping"))
    }
    return(c(none, regime = "no tipping (gradual response)"))
  }
  eta_crit <- refine(a, TRUE)
  if (eta_crit <= 0) return(c(none, regime = "cod small even without fishing"))
  # The state reached after tipping: largest stable equilibrium below the fold
  after <- equilibria(eta_crit + 1e-3 * max(1, eta_crit), p, dynamic) |>
    filter(stable, mu < g[a])
  if (nrow(after) == 0) {
    return(list(eta_crit = eta_crit, eta_rec = NA_real_, mu0 = mu0,
                regime = "complex (check bifurcation diagram)"))
  }
  i1 <- which.min(abs(g - max(after$mu)))
  if (!d$stable[i1]) {
    return(list(eta_crit = eta_crit, eta_rec = NA_real_, mu0 = mu0,
                regime = "complex (check bifurcation diagram)"))
  }
  r <- run(i1)
  k <- r[which.min(h[r])]
  eta_rec <- if (k > 1 && k < nrow(d) && !is.na(h[k - 1]) && !is.na(h[k + 1]))
    refine(k, FALSE) else h[k]
  regime <- if (eta_rec < 0) "irreversible trap" else
    "trap persists only while fishing continues"
  list(eta_crit = eta_crit, eta_rec = eta_rec, mu0 = mu0, regime = regime)
}

# Bifurcation diagram: the equilibrium curve (H(mu), mu), as line segments
# between neighbouring grid points clipped to 0 <= eta <= eta_max (so that
# branches along which mu hardly changes with eta are still drawn). Each
# segment is labelled as stable, unstable, or cod extinct.
bif_data <- function(p, dynamic, eta_max) {
  g <- mu_grid(p)
  n <- length(g)
  # Clip segment i (from point i to i + 1) of the curve (eta, mu)
  clip <- function(eta, ok) {
    e1 <- eta[-n]
    e2 <- eta[-1]
    ok <- ok[-n] & ok[-1]
    de <- e2 - e1
    t0 <- pmax(0, ifelse(de > 0, -e1 / de, (eta_max - e1) / de))
    t1 <- pmin(1, ifelse(de > 0, (eta_max - e1) / de, -e1 / de))
    flat <- !is.na(de) & de == 0
    t0[flat] <- 0
    t1[flat] <- ifelse(e1[flat] >= 0 & e1[flat] <= eta_max, 1, -1)
    i <- which(ok & !is.na(t0) & !is.na(t1) & t0 <= t1)
    tibble(i = i,
           x1 = e1[i] + t0[i] * de[i], x2 = e1[i] + t1[i] * de[i],
           y1 = g[i] + t0[i] * (g[i + 1] - g[i]),
           y2 = g[i] + t1[i] * (g[i + 1] - g[i]),
           slope = de[i])
  }
  # Branch with cod present
  B <- dmu_eq(g, 1, p, dynamic, clamp = FALSE) -
    dmu_eq(g, 0, p, dynamic, clamp = FALSE)
  h <- H_fn(g, p, dynamic)
  ok <- is.finite(h) & c(sign(B[-n]) == sign(B[-1]), TRUE) &
    c(TRUE, sign(B[-n]) == sign(B[-1]))
  s <- clip(h, ok)
  if (nrow(s) > 0) {
    # Midpoint of each segment, and the densities there
    xm <- (s$x1 + s$x2) / 2
    ym <- (s$y1 + s$y2) / 2
    N <- eco_eq(ym, xm, p, dynamic)
    s$status <- ifelse(s$slope * B[s$i] > 0, "stable", "unstable")
    extinct <- is.na(N$NC) | N$NC <= 0
    s$status[extinct] <- "cod extinct"
    if (dynamic) {
      # Stability from the full Jacobian; where a density would be negative,
      # the coexistence solution is not an equilibrium at all
      both <- !extinct & N$NF > 0
      s$status[both] <- ifelse(max_eigen(ym[both], N$NC[both], N$NF[both],
                                         xm[both], p) < 0, "stable", "unstable")
      s <- s[both, ]
    }
    s$branch <- "with cod"
  }
  if (dynamic) { # Branch without cod (flounder alone), where cod cannot invade
    h0 <- H_absent(g, p)
    s0 <- clip(h0, is.finite(h0))
    if (nrow(s0) > 0) {
      xm <- (s0$x1 + s0$x2) / 2
      ym <- (s0$y1 + s0$y2) / 2
      vC <- p$sC^2
      no_invasion <- growth_term(ym, vC, p) - fishing_term(ym, vC, xm, p) -
        hypoxia_term(ym, vC, p) -
        alpha_fn(ym, 0, vC, p$sF^2, p$alphaI, p) * flounder_alone(xm, p) <= 0
      s <- bind_rows(s, s0[no_invasion, ] |>
                       mutate(status = "cod extinct", branch = "without cod"))
    }
  }
  if (nrow(s) == 0) {
    return(tibble(grp = integer(), status = factor(), eta = numeric(),
                  mu = numeric()))
  }
  # Join consecutive segments into paths (so that dashed lines show as such)
  s |>
    arrange(branch, i) |>
    mutate(grp = cumsum(coalesce(i != lag(i) + 1 | status != lag(status) |
                                   branch != lag(branch), TRUE))) |>
    group_by(grp, status) |>
    reframe(eta = c(x1, last(x2)), mu = c(y1, last(y2))) |>
    mutate(status = factor(status, levels = c("stable", "unstable",
                                              "cod extinct")))
}


# Plot helpers ---------------------------------------------------------------

# Graphical parameters shared by the plots
base_par <- function(...) {
  par(mar = c(4.2, 4.6, 2, 1), mgp = c(2.8, 0.7, 0), las = 1, cex = 1.15,
      tcl = -0.3, ...)
}

# Reference lines: optimal size z* (dashed) and flounder mean size (dotted),
# labelled at the top of the plot (each label on the side of its line facing
# away from the other line, so that the two never overlap)
mark_sizes <- function(p) {
  abline(v = p$zstar, lty = "dashed", col = "grey50")
  abline(v = 0, lty = "dotted", col = col_flounder, lwd = 2)
  right <- p$zstar >= 0 # Is z* to the right of the flounder mean?
  top <- par("usr")[4]
  text(0, top, "flounder", adj = c(if (right) 1.08 else -0.08, 1.4),
       cex = 0.95, col = "darkseagreen4")
  text(p$zstar, top, expression(italic(z) * "*"),
       adj = c(if (right) -0.25 else 1.25, 1.4), cex = 1.05, col = "grey30")
}

# Axis range of the values y, with room at the top for the labels above
ylim_room <- function(y) {
  r <- range(y, finite = TRUE)
  r + c(0, 0.12 * diff(r))
}

# Set up a layout of plot panels (in a row) with a strip for the legend below
layout_with_legend <- function(n_panels = 1) {
  layout(matrix(c(seq_len(n_panels), rep(n_panels + 1, n_panels)), 2,
                byrow = TRUE), heights = c(1, lcm(1.5)))
}

# The legend, drawn in one row in its own strip below the panels. If it is too
# wide, the spacing is tightened and the text shrunk just enough to fit.
legend_strip <- function(...) {
  par(mar = c(0, 0, 0, 0))
  plot.new()
  avail <- 0.98 * diff(par("usr")[1:2])
  draw <- function(cex, compact = FALSE, plot = TRUE) {
    legend("center", bty = "n", horiz = TRUE, cex = cex, plot = plot,
           x.intersp = if (compact) 0.5 else 1,
           seg.len = if (compact) 1.5 else 2, ...)
  }
  width <- function(cex, compact) draw(cex, compact, plot = FALSE)$rect$w
  if (width(1, FALSE) <= avail) return(invisible(draw(1)))
  # Legend width is proportional to the text size; refine once for rounding
  cex <- avail / width(1, TRUE)
  cex <- cex * avail / width(cex, TRUE)
  invisible(draw(min(1, max(cex, 0.5)), compact = TRUE))
}


# UI -------------------------------------------------------------------------

param_slider <- function(id, label, min, max, step) {
  sliderInput(id, HTML(label), min = min, max = max, value = defaults[[id]],
              step = step, width = "100%")
}

sidebar_ui <- sidebar(
  width = 340,
  radioButtons("model", "Model",
               c("Fixed flounder density" = "fixed",
                 "Dynamic flounder density" = "dynamic")),
  actionButton("reset", "Reset to default parameters", class = "btn-sm"),
  accordion(
    open = c("Fishing", "Simulation"),
    accordion_panel(
      "Fishing",
      param_slider("eta", "Fishing effort on cod, \u03b7", 0, 30, 0.05),
      conditionalPanel(
        "input.model == 'dynamic'",
        param_slider("etaF", "Fishing effort on flounder, \u03b7<sub>F</sub>", 0, 30,
                     0.05)
      ),
      param_slider("phi", "Fishing inflection point, \u03c6", -2, 3, 0.05),
      param_slider("tau", "Fishing transition width, \u03c4", 0.1, 4, 0.05)
    ),
    accordion_panel(
      "Simulation",
      param_slider("tmax", "Simulation length", 10, 1000, 10),
      param_slider("fishing", "Fishing period", 0, defaults$tmax, 1),
      checkboxInput("burnin", "Start from the equilibrium without fishing",
                    value = defaults$burnin),
      param_slider("muC0", "Initial cod mean size, \u03bc<sub>C</sub>(0)", -3, 3, 0.05),
      param_slider("NC0", "Initial cod density, <i>N</i><sub>C</sub>(0)", 0.1, 30, 0.1)
    ),
    accordion_panel(
      "Flounder",
      param_slider("NF", "Flounder density, <i>N</i><sub>F</sub>", 0, 30, 0.1)
    ),
    accordion_panel(
      "Intrinsic growth",
      param_slider("rho", "Maximum intrinsic growth, \u03f1", 0.5, 20, 0.1),
      param_slider("theta", "Growth range width, \u03b8", 1, 10, 0.1),
      param_slider("zstar", "Optimal body size, <i>z</i>*", -2, 3, 0.05)
    ),
    accordion_panel(
      "Hypoxia",
      param_slider("kappa", "Hypoxic effect intensity, \u03ba", 0, 10, 0.05),
      param_slider("zeta", "Hypoxia inflection point, \u03b6", -2, 3, 0.05),
      param_slider("nu", "Hypoxia transition width, \u03bd", 0.1, 4, 0.05)
    ),
    accordion_panel(
      "Competition",
      param_slider("alpha0", "Baseline competition, \u03b1<sub>0</sub>", 0, 2, 0.001),
      param_slider("alphaI", "Non-size niche overlap, \u03b1<sub>I</sub>", 0, 1, 0.01),
      param_slider("w", "Competition width, \u03c9", 0.05, 3, 0.01)
    ),
    accordion_panel(
      "Trait variation",
      param_slider("sC", "Cod size SD, \u03c3<sub>C</sub>", 0.05, 2, 0.01),
      param_slider("sF", "Flounder size SD, \u03c3<sub>F</sub>", 0.05, 2, 0.01),
      param_slider("h2", "Cod heritability, <i>h</i><sup>2</sup><sub>C</sub>", 0.01, 1, 0.01)
    )
  )
)

landscape_ui <- nav_panel(
  "Model overview",
  uiOutput("landscape_note"),
  layout_columns(
    col_widths = c(6, 6), fill = FALSE,
    card(card_header(HTML(paste0(
           "Rate of change of cod mean size ",
           "<span class='text-muted fw-normal'>(&#9679; stable, ",
           "&#9675; unstable equilibrium)</span>"))),
         plotOutput("dmu_plot", height = 330)),
    card(card_header("Bifurcation diagram"),
         plotOutput("bif_plot", height = 330),
         param_slider("eta_max_bif", "Largest fishing effort shown", 1, 50, 1))
  ),
  card(card_header("Terms of the model (cod)"),
       plotOutput("terms_plot", height = 400), fill = FALSE)
)

# The animation is drawn in the browser by www/anim.js, so that playing,
# pausing and scrubbing never wait for the server.
simulation_ui <- nav_panel(
  "Simulation",
  tags$head(
    tags$script(src = "anim.js"),
    tags$style(HTML("
      #anim-title { font-size: 1.1rem; }
      #anim-gauss { width: 100%; height: 320px; display: block; }
      .ts-wrap { position: relative; }
      .ts-cursor { position: absolute; width: 0; display: none;
                   border-left: 1.5px solid #4d4d4d; pointer-events: none; }
    "))
  ),
  card(
    fill = FALSE,
    div(id = "anim-failed", style = "display: none;",
        p("The simulation failed for these parameters.")),
    div(
      id = "anim-display",
      div(class = "d-flex align-items-center gap-3",
          tags$button(id = "anim-play", type = "button",
                      class = "btn btn-primary", icon("play"), "Play"),
          tags$input(id = "anim-time", type = "range", class = "form-range",
                     min = 0, max = 1000, step = 1, value = 0,
                     style = "flex: 1;"),
          tags$select(id = "anim-speed", class = "form-select",
                      style = "width: auto;",
                      tags$option(value = "0.25", "0.25\u00d7"),
                      tags$option(value = "0.5", "0.5\u00d7"),
                      tags$option(value = "1", selected = NA, "1\u00d7"),
                      tags$option(value = "2", "2\u00d7"),
                      tags$option(value = "4", "4\u00d7"))),
      div(id = "anim-title", class = "mt-3"),
      div(id = "anim-subtitle", class = "text-muted"),
      tags$svg(id = "anim-gauss")
    )
  ),
  card(div(class = "ts-wrap",
           plotOutput("ts_plot", height = 440),
           div(id = "ts-cursor-1", class = "ts-cursor"),
           div(id = "ts-cursor-2", class = "ts-cursor")),
       fill = FALSE)
)

# Documentation, written in about.md (in the app directory)
about_ui <- nav_panel(
  "About",
  tags$style(HTML("
    .about-text { max-width: 60rem; padding: 0.5rem 1rem 2rem; }
    .about-text h2 { font-size: 1.5rem; margin-bottom: 1rem; }
    .about-text h3 { font-size: 1.2rem; margin-top: 2rem; }
    .about-text table { margin: 1rem 0; }
    .about-text th, .about-text td { padding: 0.3rem 0.8rem;
                                     border-bottom: 1px solid #dee2e6;
                                     vertical-align: top; }
  ")),
  div(class = "about-text", includeMarkdown("about.md"))
)

ui <- page_sidebar(
  title = "Cod-flounder eco-evolutionary model",
  sidebar = sidebar_ui,
  navset_card_underline(landscape_ui, simulation_ui, about_ui, id = "tabs")
)


# Server ---------------------------------------------------------------------

server <- function(input, output, session) {

  par_ids <- c("rho", "theta", "zstar", "eta", "etaF", "phi", "tau", "kappa",
               "zeta", "nu", "alpha0", "alphaI", "w", "sC", "sF", "h2", "NF")

  params <- reactive({
    p <- map(set_names(par_ids), \(id) input[[id]])
    if (input$model == "fixed") p$etaF <- 0
    p
  }) |> debounce(300)

  dynamic <- reactive(input$model == "dynamic")

  # Keep dependent controls in sync
  observeEvent(input$tmax, {
    f <- pmin(input$fishing, input$tmax)
    updateSliderInput(session, "fishing", max = input$tmax, value = f)
  })
  observeEvent(input$model, {
    updateSliderInput(session, "NF", label = if (dynamic())
      HTML("Initial flounder density, <i>N</i><sub>F</sub>(0)") else
      HTML("Flounder density, <i>N</i><sub>F</sub>"))
  })
  observeEvent(input$reset, {
    for (id in names(defaults)) {
      if (id == "burnin") updateCheckboxInput(session, id, value = defaults[[id]])
      else updateSliderInput(session, id, value = defaults[[id]])
    }
    updateSliderInput(session, "fishing", max = defaults$tmax)
  })

  # Analysis of the model version shown
  folds <- reactive(fold_points(params(), dynamic()))

  output$landscape_note <- renderUI({
    if (!dynamic()) return(NULL)
    p(class = "text-muted", HTML(paste0(
      "Dynamic flounder density: everything on this tab assumes that cod and ",
      "flounder densities are at their ecological equilibrium for the given ",
      "cod mean size (so the initial flounder density plays no role here). ",
      "Flounder are fished in proportion to cod, at the ratio ",
      "&eta;<sub>F</sub>/&eta; set in the sidebar.")))
  })

  output$dmu_plot <- renderPlot({
    p <- params()
    dyn <- dynamic()
    g <- seq(min(-3, p$zstar - 2), max(3, p$zstar + 2), length.out = 801)
    curves <- list(dmu_eq(g, 0, p, dyn), dmu_eq(g, p$eta, p, dyn))
    eqs <- list(equilibria(0, p, dyn), equilibria(p$eta, p, dyn))
    cols <- c("grey40", col_fishing)
    layout_with_legend()
    base_par()
    plot(NA, xlim = range(g), ylim = ylim_room(unlist(curves)), xaxs = "i",
         xlab = expression("Cod mean body size," ~ mu[C]),
         ylab = expression(d * mu[C] / d * t))
    abline(h = 0, col = "grey60")
    mark_sizes(p)
    for (k in 1:2) {
      lines(g, curves[[k]], col = cols[k], lwd = 2.2)
      points(eqs[[k]]$mu, rep(0, nrow(eqs[[k]])), col = cols[k], bg = "white",
             pch = ifelse(eqs[[k]]$stable, 16, 21), cex = 1.4, lwd = 1.6)
    }
    legend_strip(legend = c("no fishing", sprintf("fishing, \u03b7 = %.2f", p$eta)),
                 col = cols, lwd = 2.2)
  }, execOnResize = TRUE) # Re-run on resize, as legends etc. are fitted

  output$bif_plot <- renderPlot({
    p <- params()
    eta_max <- input$eta_max_bif
    d <- bif_data(p, dynamic(), eta_max)
    validate(need(nrow(d) > 0, "No equilibria in this range of fishing effort."))
    f <- folds()
    styles <- list(stable = list(col = col_cod, lty = "solid"),
                   unstable = list(col = col_cod, lty = "dashed"),
                   `cod extinct` = list(col = "grey45", lty = "dotted"))
    layout_with_legend()
    base_par()
    plot(NA, xlim = c(0, eta_max), ylim = range(d$mu, finite = TRUE), xaxs = "i",
         xlab = expression("Fishing effort," ~ eta),
         ylab = expression("Equilibrium" ~ mu[C]))
    abline(h = 0, lty = "dotted", col = col_flounder, lwd = 2)
    abline(v = p$eta, col = col_fishing, lwd = 1.5)
    # Tipping points
    tips <- list(list(eta = f$eta_crit, label = expression(eta[crit])),
                 list(eta = f$eta_rec, label = expression(eta[rec])))
    for (tp in tips) {
      if (is.na(tp$eta) || tp$eta < 0 || tp$eta > eta_max) next
      abline(v = tp$eta, lty = "dashed", col = "grey40")
      text(tp$eta, par("usr")[4], tp$label, adj = c(-0.15, 1.3), col = "grey30",
           cex = 1.15)
    }
    for (s in split(d, d$grp)) {
      st <- styles[[as.character(s$status[1])]]
      lines(s$eta, s$mu, col = st$col, lty = st$lty, lwd = 2.2)
    }
    shown <- names(styles)[names(styles) %in% d$status]
    legend_strip(legend = shown, lwd = 2.2,
                 col = map_chr(styles[shown], "col"),
                 lty = map_chr(styles[shown], "lty"))
  }, execOnResize = TRUE) # Re-run on resize, as legends etc. are fitted

  # The terms of the model, each with the sign it enters with, superimposed
  # with their sum: the per capita growth rate of cod when rare (i.e. without
  # competition among cod), and the selection on cod mean size, whose sum is
  # d mu_C / dt divided by h2 (so its zeros are the equilibria shown above)
  output$terms_plot <- renderPlot({
    p <- params()
    vC <- p$sC^2
    vF <- p$sF^2
    m <- seq(min(-3, p$zstar - 2), max(3, p$zstar + 2), length.out = 401)
    # Flounder density: fixed, or at ecological equilibrium (current eta)
    NF <- eco_eq(m, p$eta, p, dynamic(), clamp = TRUE)$NF
    growth <- list(growth = growth_term(m, vC, p),
                   fishing = -fishing_term(m, vC, p$eta, p),
                   hypoxia = -hypoxia_term(m, vC, p),
                   competition = -alpha_fn(m, 0, vC, vF, p$alphaI, p) * NF)
    selection <- list(growth = growth_sel(m, vC, p),
                      fishing = fishing_sel(m, vC, p$eta, p),
                      hypoxia = hypoxia_sel(m, vC, p),
                      competition = comp_sel(m, 0, vC, vF, NF, p$alphaI, p))
    cols <- c(growth = "steelblue", fishing = col_fishing,
              hypoxia = "darkorange2", competition = "seagreen4",
              total = "black")
    lwds <- c(growth = 2, fishing = 2, hypoxia = 2, competition = 2, total = 3.2)
    panel <- function(terms, title) {
      curves <- c(terms, list(total = Reduce(`+`, terms)))
      plot(NA, xlim = range(m), ylim = ylim_room(unlist(curves)), xaxs = "i",
           xlab = expression("Cod mean body size," ~ mu[C]), ylab = "")
      abline(h = 0, col = "grey60")
      mark_sizes(p)
      for (k in names(curves)) lines(m, curves[[k]], col = cols[[k]], lwd = lwds[[k]])
      title(main = title, adj = 0, font.main = 1, cex.main = 1)
    }
    layout_with_legend(2)
    base_par()
    panel(growth, "Per capita growth rate of rare cod")
    panel(selection, expression("Selection on body size," ~
                                  (d * mu[C] / d * t) / italic(h)[C]^2))
    legend_strip(legend = c("intrinsic growth",
                            sprintf("fishing (\u03b7 = %.2f)", p$eta),
                            "hypoxia", "competition with flounder", "total"),
                 col = cols, lwd = lwds)
  }, execOnResize = TRUE) # Re-run on resize, as legends etc. are fitted

  # Simulation
  sim <- reactive({
    p <- params()
    y0 <- c(NC = input$NC0, muC = input$muC0)
    if (dynamic()) y0 <- c(y0, NF = p$NF)
    tryCatch(
      simulate(p, dynamic(), input$fishing[1], input$fishing[2], input$tmax,
               y0, input$burnin),
      error = \(e) NULL)
  })

  # Animation: send the trajectory to the browser (www/anim.js), which takes
  # care of playback. Only done while the Simulation tab is open, so that
  # parameter changes on the other tab do not trigger simulations.
  observe({
    req(input$tabs == "Simulation")
    s <- sim()
    p <- params()
    if (is.null(s)) {
      session$sendCustomMessage("anim-data", list(failed = TRUE))
      return()
    }
    spread <- 3 * max(p$sC, p$sF)
    session$sendCustomMessage("anim-data", list(
      time = s$time, NC = s$NC, muC = s$muC, NF = s$NF,
      tmax = max(s$time), fishing = input$fishing, eta = p$eta,
      sC = p$sC, sF = p$sF, zstar = p$zstar,
      zlim = c(min(s$muC, 0, p$zstar) - spread, max(s$muC, 0, p$zstar) + spread),
      ymax = max(s$NC / (p$sC * sqrt(2*pi)), s$NF / (p$sF * sqrt(2*pi)), 1e-6),
      col_cod = col_cod, col_flounder = col_flounder, col_fishing = col_fishing
    ))
  })

  # Time series in base R (fast). The browser draws the moving time cursor on
  # top, so it needs to know where the panels sit on the image (as fractions
  # of the image size, measured from the bottom left).
  output$ts_plot <- renderPlot({
    s <- sim()
    validate(need(!is.null(s), "The simulation failed for these parameters."))
    xlim <- range(s$time)
    fishing_band <- function() {
      u <- par("usr")
      rect(input$fishing[1], u[3], input$fishing[2], u[4], border = NA,
           col = adjustcolor(col_fishing, alpha.f = 0.08))
    }
    par(mfrow = c(2, 1), mar = c(2, 5, 2, 1), oma = c(2.5, 0, 0, 0),
        mgp = c(3.4, 0.7, 0), las = 1, cex = 1.25, tcl = -0.3)
    plot(NA, xlim = xlim, ylim = range(s$muC, 0), xaxs = "i",
         xlab = "", ylab = expression("Mean size," ~ mu))
    fishing_band()
    abline(h = 0, col = col_flounder, lwd = 2.5)
    lines(s$time, s$muC, col = col_cod, lwd = 2.5)
    title(main = "Mean body size", adj = 0, font.main = 1, cex.main = 1)
    legend("topright", c("cod", "flounder"), col = c(col_cod, col_flounder),
           lwd = 2.5, horiz = TRUE, bty = "n", cex = 0.9)
    y1 <- grconvertY(par("usr")[3:4], "user", "ndc")
    x <- grconvertX(xlim, "user", "ndc")
    # Densities, with some room at the top for labelling the fishing period
    dens <- cbind(s$NC, s$NF)
    plot(NA, xlim = xlim, ylim = c(0, 1.18 * max(dens)), xaxs = "i",
         xlab = "", ylab = "Density")
    fishing_band()
    lines(s$time, dens[, 2], col = col_flounder, lwd = 2.5)
    lines(s$time, dens[, 1], col = col_cod, lwd = 2.5)
    f <- pmin(pmax(input$fishing, xlim[1]), xlim[2])
    labels <- tibble(from = c(xlim[1], f[1], f[2]), to = c(f[1], f[2], xlim[2]),
                     on = c(FALSE, TRUE, FALSE)) |>
      mutate(text = if_else(on, "Fishing on", "Fishing off"))
    labels <- labels[1.15 * strwidth(labels$text, cex = 0.85) <
                       labels$to - labels$from, ]
    # "Fishing on" centred in its region; "Fishing off" right where its
    # region begins
    on <- labels[labels$on, ]
    off <- labels[!labels$on, ]
    text((on$from + on$to) / 2, par("usr")[4], on$text, adj = c(0.5, 1.5),
         cex = 0.85, col = col_fishing)
    text(off$from + 0.01 * diff(xlim), par("usr")[4], off$text,
         adj = c(0, 1.5), cex = 0.85, col = "grey40")
    title(main = "Population density", adj = 0, font.main = 1, cex.main = 1)
    mtext("Time", side = 1, outer = TRUE, line = 1.2, cex = 1.3)
    y2 <- grconvertY(par("usr")[3:4], "user", "ndc")
    session$sendCustomMessage("ts-map", list(t = xlim, x = x, y1 = y1, y2 = y2))
  }, execOnResize = TRUE) # Re-run on resize, as legends etc. are fitted
}

shinyApp(ui, server)
