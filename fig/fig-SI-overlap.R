library(tidyverse)
library(mvtnorm)


# Create data frame with points on the perimeter of an ellipse
# Input:
# - m: Two-component vector, with the means in the x- and y-direction
# - P: Covariance matrix
# - p: What fraction of the distribution to show (0 < p < 1)
# - res: Resolution - i.e., how many points should the ellipse will consist of
# Output:
# - A tibble with two columns (x and y) and res rows; each row
#   contains the coordinates of one point of the ellipse
ellipseframe <- function(m = rep(0, 2), P = diag(2), p = 0.95, res = 101) {
  s <- -2*log(1 - p)
  eP <- eigen(P*s)
  sqrtP <- eP$vectors %*% diag(sqrt(eP$values)) %*% t(eP$vectors)
  phi <- seq(0, 2*pi, l = res)
  ellipse <- tibble(x = cos(phi), y = sin(phi))
  for (i in 1:nrow(ellipse)) ellipse[i,] <- as.list(sqrtP %*% unlist(ellipse[i,]))
  mutate(ellipse, x = x + m[1], y = y + m[2])
}


tibble(codMean = list(c(2, -1), c(0, -1), c(-2, -1)),
       flounderMean = list(c(0, 1), c(0, 1), c(0, 1))) |>
  mutate(codBody = map_dbl(codMean, \(x) x[1])) |>
  mutate(cod = map(codMean, \(m) ellipseframe(m))) |>
  mutate(flounder = map(flounderMean, \(m) ellipseframe(m))) |>
  select(!ends_with("Mean")) |>
  pivot_longer(cod | flounder, names_to = "species", values_to = "value") |>
  unnest(value) |>
  mutate(codBody = case_when(
    codBody == max(codBody) ~ "Cod larger than flounder",
    codBody == min(codBody) ~ "Flounder larger than cod",
    .default = "Equal body size"
  )) |>
  ggplot(aes(x = x, y = y, color = species, fill = species)) +
  geom_polygon(alpha = 0.2) +
  labs(x = "Body size", y = "Other trait") +
  scale_color_manual(name = "Species", values = c("cornflowerblue", "darkseagreen")) +
  scale_fill_manual(name = "Species", values = c("cornflowerblue", "darkseagreen")) +
  facet_grid(. ~ codBody) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank())
#ggsave("fig-SI-overlap.pdf", width = 7, height = 2)
