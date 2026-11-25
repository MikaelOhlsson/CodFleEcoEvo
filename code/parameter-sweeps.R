# Parameter sweeps of the tipping point in cod body size (reviewer comments
# E2, R1.3b, R1.6).
#
# Uses the reduced model of the main text (flounder density and mean size held
# fixed), so only cod mean body size mu changes:
#
#   dmu/dt = h2 * [A(mu) - eta * B(mu)],
#
# where B(mu) > 0 is the selection from fishing per unit fishing effort eta,
# and A(mu) collects all other selection (growth, hypoxia, competition).
# Because eta enters linearly, the equilibria for a given eta are the
# solutions of eta = H(mu) := A(mu) / B(mu). The curve (H(mu), mu) is thus the
# full bifurcation diagram, and an equilibrium is stable iff H'(mu) < 0.
# The saddle-node bifurcations (tipping points) are the local extrema of H:
#   - local maximum of H: eta_crit, above which the large-size state is lost;
#   - local minimum of H: eta_rec, below which the small-size state is lost
#     (if eta_rec < 0, the small state persists even without fishing, i.e.
#     the trap is irreversible by reducing fishing alone).
# Note that h2 cancels from H: heritability affects the speed of evolution,
# but not where the tipping points are. Likewise, alpha0, alphaI and Fpop
# only enter through their product, so sweeping alphaI covers all three.
#
# Regimes, for given parameters:
#   "irreversible trap"      bistable already at eta = 0 (eta_rec < 0)
#   "reversible hysteresis"  bistable for 0 < eta_rec < eta < eta_crit
#   "no tipping"             single state; cod size declines smoothly
#   "small without fishing"  single state, and cod are already smaller than
#                            flounder without any fishing
#
# Run from the repository root; figures are saved to figures/.

library(tidyverse)
library(patchwork)


# Model --------------------------------------------------------------------

# The competition term of markdown/cod-flounder-basic.qmd (which produces
# Fig. 1) has (s^2 + 2*sF^2 + w^2)^3 under the square root in its
# denominator, whereas Eq. 4 of the manuscript (and the full model in
# cod-flounder-feedback.qmd) has (s^2 + sF^2 + w^2)^3. With the Fig. 1
# parameters, this moves eta_crit from about 1.56 (qmd) to about 5.5 (Eq. 4).
# Set to TRUE to reproduce the qmd version.
use_qmd_denominator <- FALSE

# Fig. 1 parameter values (see markdown/cod-flounder-basic.qmd):
pars_fig1 <- list(
  mF     = 0,          # Flounder mean (scaled) body size
  rho    = 5,          # Maximum intrinsic growth rate
  s      = 0.5,        # Cod body size standard deviation
  sF     = 0.5,        # Flounder body size standard deviation
  theta  = 5,          # Trait range where population growth is possible
  zstar  = 1,          # Optimal body size for cod
  phi    = 0,          # Fishing inflection point
  tau    = 0.5,        # Fishing transition width
  kappa  = 2,          # Hypoxic effect intensity
  zeta   = 1,          # Hypoxic effect inflection point
  nu     = 1.5,        # Hypoxic effect transition width
  Fpop   = 6,          # Flounder population density
  w      = 0.33,       # Competition width of body size
  alpha0 = 1/sqrt(2),  # Baseline competition from body size
  alphaI = 1           # Competition from traits other than body size
)

# Selection from fishing per unit fishing effort, B(mu) (without the h2)
sel_fishing <- function(m, p) {
  with(p, s^2 / sqrt(pi*(2*s^2 + tau^2)) * exp(-(m - phi)^2 / (2*s^2 + tau^2)))
}

# All other selection, A(mu): growth, hypoxia, and competition (without h2)
sel_other <- function(m, p) {
  with(p, {
    S <- s^2 + sF^2 + w^2
    S_denom <- if (use_qmd_denominator) s^2 + 2*sF^2 + w^2 else S
    2*rho*s^2 / theta^2 * (zstar - m) -
      kappa*s^2 / sqrt(pi*(2*s^2 + nu^2)) * exp(-(m - zeta)^2 / (2*s^2 + nu^2)) +
      s^2*Fpop*alpha0*alphaI*(m - mF) / sqrt(2*pi*S_denom^3) *
      exp(-(m - mF)^2 / (2*S))
  })
}

# Fishing effort at which mu is an equilibrium
H <- function(m, p) sel_other(m, p) / sel_fishing(m, p)


# Tipping points -----------------------------------------------------------

mu_grid <- seq(-5, 5, by = 0.002)

# Locate the saddle-node bifurcations and classify the regime, for the
# parameter list p. Returns a one-row tibble.
fold_points <- function(p, grid = mu_grid) {
  h <- H(grid, p)
  # Cod mean size at the largest stable equilibrium without fishing
  # (stable equilibria are where H crosses 0 from above):
  i0 <- which(h[-length(h)] > 0 & h[-1] <= 0)
  mu0 <- uniroot(\(m) H(m, p), grid[c(max(i0), max(i0) + 1)])$root
  # Interior local extrema of H, refined numerically:
  dh <- sign(diff(h))
  ext <- which(diff(dh) != 0) + 1
  refine <- function(i, maximum) {
    optimize(\(m) H(m, p), grid[c(i - 1, i + 1)], maximum = maximum)$objective
  }
  type <- ifelse(dh[ext - 1] > 0, "max", "min") # Slope of H before extremum
  single <- if (mu0 < p$mF) "small without fishing" else "no tipping"
  if (length(ext) == 0) {
    return(tibble(eta_crit = NA_real_, eta_rec = NA_real_, mu0 = mu0,
                  regime = single))
  }
  if (length(ext) != 2 || !identical(type, c("min", "max"))) {
    return(tibble(eta_crit = NA_real_, eta_rec = NA_real_, mu0 = mu0,
                  regime = "other (check)"))
  }
  eta_rec <- refine(ext[1], maximum = FALSE)
  eta_crit <- refine(ext[2], maximum = TRUE)
  regime <- case_when(
    eta_crit <= 0 ~ single, # Bistability only at (unphysical) negative eta
    eta_rec < 0   ~ "irreversible trap",
    TRUE          ~ "reversible hysteresis"
  )
  if (eta_crit <= 0) eta_crit <- eta_rec <- NA_real_
  tibble(eta_crit = eta_crit, eta_rec = eta_rec, mu0 = mu0, regime = regime)
}

# Sweep over a table of parameter values. Each column of `grid` must be a
# parameter name; all other parameters are taken from `base`.
run_sweep <- function(grid, base = pars_fig1) {
  grid |>
    mutate(res = pmap(grid, \(...) fold_points(modifyList(base, list(...))))) |>
    unnest(res)
}

regime_levels <- c("irreversible trap", "reversible hysteresis",
                   "no tipping", "small without fishing", "other (check)")
regime_cols <- setNames(c("firebrick", "goldenrod", "steelblue", "grey60", "black"),
                        regime_levels)


# Sanity check at the Fig. 1 parameters ------------------------------------

print(fold_points(pars_fig1))


# Sweep 1: alphaI x kappa --------------------------------------------------

# How much fishing it takes to tip the system, as a function of the strength
# of competition (non-size niche overlap alphaI) and of hypoxia (kappa)

sweep_ak <-
  crossing(alphaI = seq(0.05, 1, by = 0.01), kappa = seq(0, 4, by = 0.05)) |>
  run_sweep() |>
  mutate(regime = factor(regime, levels = regime_levels))

font_size <- 12

theme_sweep <- theme_bw(base_size = font_size) +
  theme(panel.grid = element_blank())

mark_fig1 <- annotate("point", x = pars_fig1$alphaI, y = pars_fig1$kappa,
                      shape = 21, fill = "white", size = 2.5)

p_ak_regime <- ggplot(sweep_ak, aes(x = alphaI, y = kappa, fill = regime)) +
  geom_raster() +
  mark_fig1 +
  scale_fill_manual(values = regime_cols, breaks = regime_levels,
                    drop = TRUE) +
  scale_x_continuous(expand = expansion(0, 0)) +
  scale_y_continuous(expand = expansion(0, 0)) +
  labs(x = expression("Non-size niche overlap," ~ alpha[I]),
       y = expression("Hypoxic effect intensity," ~ kappa),
       fill = NULL) +
  theme_sweep +
  theme(legend.position = "bottom", legend.direction = "vertical")

p_ak_etacrit <- ggplot(sweep_ak, aes(x = alphaI, y = kappa)) +
  geom_raster(aes(fill = eta_crit)) +
  geom_contour(aes(z = eta_crit), colour = "white", linewidth = 0.3,
               breaks = c(0.5, 1, 2, 4, 8, 16), na.rm = TRUE) +
  # Boundary between irreversible traps and reversible hysteresis:
  geom_contour(aes(z = eta_rec), breaks = 0, colour = "black",
               linetype = "dashed", linewidth = 0.4, na.rm = TRUE) +
  mark_fig1 +
  scale_fill_viridis_c(trans = "log10", na.value = "grey85",
                       breaks = c(0.01, 0.1, 1, 10),
                       labels = c("0.01", "0.1", "1", "10")) +
  scale_x_continuous(expand = expansion(0, 0)) +
  scale_y_continuous(expand = expansion(0, 0)) +
  labs(x = expression("Non-size niche overlap," ~ alpha[I]),
       y = expression("Hypoxic effect intensity," ~ kappa),
       fill = expression(eta[crit])) +
  theme_sweep +
  theme(legend.position = "bottom",
        legend.key.width = unit(2.5, "lines"))

fig_ak <- p_ak_regime + p_ak_etacrit + plot_annotation(tag_levels = "A")

ggsave("figures/fig-sweep-alphaI-kappa.pdf", fig_ak, width = 9, height = 5.5)


# Sweep 2: alphaI x eta ----------------------------------------------------

# Which states are stable, for given niche overlap and fishing effort
# (hypoxia and all else at Fig. 1 values)

eta_max_plot <- 10

sweep_a <- tibble(alphaI = seq(0.05, 1, by = 0.005)) |> run_sweep()

regions_ae <-
  crossing(sweep_a, eta = seq(0, eta_max_plot, l = 201)) |>
  mutate(state = case_when(
    is.na(eta_crit) ~ "single state",
    eta < eta_rec   ~ "large only",
    eta > eta_crit  ~ "small only",
    TRUE            ~ "bistable"
  )) |>
  mutate(state = factor(state, levels = c("large only", "bistable",
                                          "small only", "single state")))

p_ae <- ggplot(regions_ae, aes(x = alphaI, y = eta)) +
  geom_raster(aes(fill = state)) +
  geom_line(data = sweep_a, aes(y = eta_crit), na.rm = TRUE) +
  geom_line(data = filter(sweep_a, eta_rec >= 0), aes(y = eta_rec),
            linetype = "dashed", na.rm = TRUE) +
  annotate("point", x = pars_fig1$alphaI, y = 2, shape = 21,
           fill = "white", size = 2.5) + # Fishing effort used in Fig. 1
  scale_fill_manual(values = c("large only" = "steelblue",
                               "bistable" = "goldenrod",
                               "small only" = "firebrick",
                               "single state" = "grey85")) +
  scale_x_continuous(expand = expansion(0, 0)) +
  scale_y_continuous(expand = expansion(0, 0)) +
  coord_cartesian(ylim = c(0, eta_max_plot)) +
  labs(x = expression("Non-size niche overlap," ~ alpha[I]),
       y = expression("Fishing effort," ~ eta),
       fill = NULL) +
  theme_sweep

ggsave("figures/fig-sweep-alphaI-eta.pdf", p_ae, width = 6, height = 4.5)


# Sweep 3: one parameter at a time -----------------------------------------

# Range of each remaining parameter (alphaI, kappa covered above; h2 does not
# affect the tipping points; alpha0 and Fpop act only through alphaI)
oat_ranges <- list(
  rho   = seq(2.5, 10, l = 101),
  theta = seq(3, 8, l = 101),
  zstar = seq(0.5, 2, l = 101),
  s     = seq(0.25, 1, l = 101),
  sF    = seq(0.25, 1, l = 101),
  phi   = seq(-0.5, 1.5, l = 101),
  tau   = seq(0.25, 2, l = 101),
  zeta  = seq(0, 2, l = 101),
  nu    = seq(0.5, 3, l = 101),
  w     = seq(0.1, 1, l = 101)
)

oat_labels <- c(
  rho   = "Max.~intrinsic~growth~('\u03F1')",
  theta = "Growth~range~width~(theta)",
  zstar = "Optimal~body~size~(italic(z)*'*')",
  s     = "Cod~size~SD~(sigma[C])",
  sF    = "Flounder~size~SD~(sigma[F])",
  phi   = "Fishing~inflection~point~(phi)",
  tau   = "Fishing~transition~width~(tau)",
  zeta  = "Hypoxia~inflection~point~(zeta)",
  nu    = "Hypoxia~transition~width~(nu)",
  w     = "Competition~width~(omega)"
)

sweep_oat <-
  imap(oat_ranges, \(vals, par) {
    tibble(!!par := vals) |>
      run_sweep() |>
      rename(value = all_of(par)) |>
      mutate(parameter = par)
  }) |>
  list_rbind() |>
  mutate(regime = factor(regime, levels = regime_levels),
         parameter = factor(parameter, levels = names(oat_ranges),
                            labels = oat_labels[names(oat_ranges)]))

defaults_oat <-
  tibble(parameter = factor(oat_labels, levels = oat_labels),
         value = unlist(pars_fig1[names(oat_ranges)]))

p_oat <- ggplot(sweep_oat, aes(x = value)) +
  geom_ribbon(aes(ymin = pmax(eta_rec, 0), ymax = eta_crit), na.rm = TRUE,
              fill = "goldenrod", alpha = 0.4) +
  geom_line(aes(y = eta_crit), na.rm = TRUE) +
  geom_line(aes(y = pmax(eta_rec, 0)), linetype = "dashed", na.rm = TRUE) +
  geom_rug(aes(colour = regime), sides = "b", length = unit(0.05, "npc")) +
  geom_vline(data = defaults_oat, aes(xintercept = value),
             linetype = "dotted", colour = "grey40") +
  facet_wrap(~ parameter, scales = "free", ncol = 4,
             labeller = label_parsed, strip.position = "bottom") +
  scale_colour_manual(values = regime_cols, breaks = regime_levels,
                      drop = TRUE) +
  labs(x = NULL, y = expression("Fishing effort," ~ eta), colour = NULL) +
  theme_sweep +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "bottom")

ggsave("figures/fig-sweep-oat.pdf", p_oat, width = 10, height = 7,
       device = cairo_pdf)

