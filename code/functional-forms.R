# Schematic figure of the functional forms in the eco-evolutionary model.
#
# Top row: the terms under the underbraces of Eq. 3, as functions of the mean
# trait: intrinsic growth (parabola), fishing / hypoxic mortality (sigmoid),
# and resource competition (Gaussian).
# Bottom row: the corresponding terms under the underbraces of Eq. 4. Each
# equals sigma^2 times the derivative of the term directly above it (up to the
# sign with which the term enters the equations).
# In each column, the width parameter is varied (colours); all other parameters
# are set to their Fig. 1 values. Dashed vertical lines mark the location
# parameters, dotted horizontal lines the height parameters.
#
# Run from the repository root; the figure is saved to figures/.

library(tidyverse)
library(patchwork)


Q <- function(z) pnorm(sqrt(2) * z) # (erf(z) + 1) / 2, as in the manuscript

# Fig. 1 parameter values (not varied here):
v <- 0.5^2 # Trait variance sigma^2
rho <- 5 # Maximum intrinsic growth rate
zstar <- 1 # Optimal body size
eta <- 1 # Fishing effort (or kappa: maximum hypoxic mortality)
phi <- 1.1 # Midpoint of fishing (or zeta: of hypoxic) mortality

cols <- c("steelblue", "goldenrod", "forestgreen")
mu_range <- seq(-4, 6, l = 401)
dmu_range <- seq(-5, 5, l = 401)
font_size <- 12 # Base font size (pt), also used for in-panel annotations


# Data ---------------------------------------------------------------------

# Intrinsic growth, varying theta
d_growth <- crossing(mu = mu_range, par = c(3, 5, 8)) |>
  mutate(term = rho * (par^2 - (zstar - mu)^2 - v) / par^2,
         sel  = 2 * rho * v * (zstar - mu) / par^2)

# Fishing or hypoxic mortality, varying tau (or nu)
d_mort <- crossing(mu = mu_range, par = c(0.5, 2, 4)) |>
  mutate(s    = sqrt(2 * v + par^2),
         term = eta * Q((mu - phi) / s),
         sel  = eta * v * exp(-(mu - phi)^2 / s^2) / sqrt(pi * s^2))

# Resource competition per unit alpha_0 N_j (with alpha_I = 1), varying omega;
# the competitor's trait variance is set equal to the focal species'
d_comp <- crossing(mu = dmu_range, par = c(0.5, 1, 2)) |>
  mutate(s2   = 2 * v + par^2,
         term = exp(-mu^2 / (2 * s2)) / sqrt(2 * pi * s2),
         sel  = mu * v / s2 * term)


# Plotting -----------------------------------------------------------------

plot_fn <- function(d, y, xlab, ylab, collab, legend = TRUE) {
  ggplot(d, aes(x = mu, y = {{y}}, colour = as_factor(par))) +
    geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.3) +
    geom_line(linewidth = 0.7) +
    scale_colour_manual(values = cols) +
    labs(x = xlab, y = ylab, colour = collab) +
    theme_bw(base_size = font_size) +
    theme(panel.grid = element_blank(),
          legend.position = if (legend) "top" else "none",
          legend.margin = margin(0, 0, -4, 0),
          legend.key.width = unit(1.2, "lines"))
}

# Dashed vertical line and label marking a location parameter
mark_location <- function(x, label, y, hjust = -0.15) {
  list(geom_vline(xintercept = x, linetype = "dashed", colour = "grey50"),
       annotate("text", x = x, y = y, label = label, parse = TRUE,
                hjust = hjust, colour = "grey30", size = font_size / .pt))
}

# Dotted horizontal line and label marking a height parameter
mark_height <- function(y, label, x, vjust = 1.4, parse = TRUE) {
  list(geom_hline(yintercept = y, linetype = "dotted", colour = "grey50"),
       annotate("text", x = x, y = y, label = label, parse = parse,
                hjust = 0, vjust = vjust, colour = "grey30",
                size = font_size / .pt))
}

mu_lab <- expression("Mean body size," ~ mu[i])
dmu_lab <- expression("Body size difference," ~ mu[i] - mu[j])
mort_par <- expression("Width," ~ tau ~ "(or" ~ nu * ")")

p_growth <- d_growth |>
  plot_fn(term, mu_lab, "Intrinsic growth", expression("Width," ~ theta)) +
  mark_location(zstar, "italic(z) * \"*\"", y = -7.5) +
  # Plotmath has no varrho, so use the Unicode character:
  mark_height(rho, "\u03F1", x = min(mu_range), parse = FALSE)

p_mort <- d_mort |>
  plot_fn(term, mu_lab, "Fishing (or hypoxic) mortality", mort_par) +
  mark_location(phi, "phi ~ \"(or\" ~ zeta * \")\"", y = 0.08) +
  mark_height(eta, "eta ~ \"(or\" ~ kappa * \")\"", x = min(mu_range))

p_comp <- d_comp |>
  plot_fn(term, dmu_lab, expression("Resource competition," ~
                                      alpha[ij] / alpha[0]),
          expression("Width," ~ omega))

p_growth_sel <- d_growth |>
  plot_fn(sel, mu_lab, "Selection from growth", NULL, legend = FALSE) +
  mark_location(zstar, "italic(z) * \"*\"", y = -1.2)

p_mort_sel <- d_mort |>
  plot_fn(sel, mu_lab, "Selection from fishing (or hypoxia)", NULL,
          legend = FALSE) +
  mark_location(phi, "phi ~ \"(or\" ~ zeta * \")\"", y = 0.14, hjust = -0.6)

p_comp_sel <- d_comp |>
  plot_fn(sel, dmu_lab, "Selection from competition", NULL, legend = FALSE)


fig <- (p_growth + p_mort + p_comp) /
  (p_growth_sel + p_mort_sel + p_comp_sel) +
  plot_annotation(tag_levels = "A")

# cairo_pdf is needed to render the Unicode varrho:
ggsave("figures/fig-functional-forms.pdf", fig, width = 9.5, height = 6,
       device = cairo_pdf)
