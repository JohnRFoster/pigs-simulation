library(tidyverse)
library(ggeffects)
library(effects)
library(marginaleffects)
library(mgcv)

path <- "analysis/GAMs"
y <- "mpe"

model_rds <- paste0(y, "-allSingle-null.rds")
rds <- read_rds(file.path(path, model_rds))
m0 <- rds$fit

model_rds <- paste0(y, "-allSingle.rds")
rds <- read_rds(file.path(path, model_rds))
m1 <- rds$fit

anova.gam(m0, m1, test = "Chisq")

model_rds <- paste0(y, "-allInteractions.rds")
rds <- read_rds(file.path(path, model_rds))
m2 <- rds$fit

anova.gam(m1, m2, test = "Chisq")

model_rds <- paste0(y, "-cut1.rds")
rds <- read_rds(file.path(path, model_rds))
m3 <- rds$fit
data <- rds$data

anova.gam(m2, m3, test = "Chisq")

summary(m3)
fit <- m3

terms <- attr(fit$terms, "term.labels")
labels <- gsub("_", " ", terms)
labels <- gsub("I ", "", labels)

ylab <- if_else(y == "bias", y, paste0("log(", y, ")"))

i <- 2
x <- terms[i]

my_theme <- function(s = 12) {
  theme_bw() +
    theme(
      axis.title = element_text(size = s),
      axis.text = element_text(size = s - 4),
      title = element_text(size = s + 2)
    )
}

save_gg <- function(dest, gg, path) {
  if (!dir.exists(path)) {
    dir.create(path, showWarnings = FALSE, recursive = TRUE)
  }
  ggsave(
    filename = dest,
    plot = gg,
    device = "jpeg",
    path = path,
    width = 7,
    height = 5,
    units = "in",
    dpi = "retina",
    bg = "white"
  )
}

plot_marginal_effects <- function(term, label, ylab) {
  g <- plot_predictions(
    model = fit,
    condition = term,
    rug = TRUE,
    draw = TRUE
  ) +
    labs(
      title = paste0("Predictions for different values of ", label),
      subtitle = "All other covariates held at their means",
      x = paste0(label, " (standardized)"),
      y = ylab
    ) +
    my_theme()

  dest <- paste0(ylab, "_marginalPrediction_", term, ".jpeg")
  save_gg(dest, g, "Plots/GAMs")

  g <- plot_slopes(fit, variables = term, condition = term) +
    labs(
      title = paste0(
        "Marginal effect of ",
        label,
        " for different values of ",
        label
      ),
      x = paste0(label, " (standardized)")
    ) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    my_theme()

  dest <- paste0(ylab, "_marginalSlopes_", term, ".jpeg")
  save_gg(dest, g, "Plots/GAMs")
}

for (i in 2:length(terms)) {
  plot_marginal_effects(terms[i], labels[i], ylab)
}


## plot with interactions using quartiles
plot_interactions <- function(term, label, ylab) {
  my_cond <- vector(mode = "list", length = 2)
  names(my_cond) <- c("property_area", term)
  my_cond[[2]] <- "quartile"
  g <- plot_predictions(
    model = fit,
    condition = my_cond
  ) +
    labs(
      title = paste0(
        "Interaction between quartiles of ",
        label,
        " across property area"
      ),
      x = " property area (standardized)",
      y = ylab,
      color = label,
      fill = label
    ) +
    my_theme(10) +
    theme(legend.position = "bottom")

  dest <- paste0(ylab, "_interactionPrediction_propertyArea_", term, ".jpeg")
  save_gg(dest, g, "Plots/GAMs")
}

for (i in 3:length(terms)) {
  plot_interactions(terms[i], labels[i], ylab)
}

pred <- predict.gam(fit, se.fit = TRUE)

fit_ci <- pred %>%
  map_dfc(as.numeric) %>%
  mutate(ymin = fit - (1.96 * se.fit), ymax = fit + (1.96 * se.fit), y = data$y)

r2 <- function(observed, predicted) cor(observed, predicted)^2
rsq <- round(r2(fit_ci$y, fit_ci$fit), 2)


g <- fit_ci |>
  ggplot() +
  aes(x = y, y = fit) +
  geom_point(size = 0.5) +
  geom_abline(aes(intercept = 0, slope = 1, color = "1:1")) +
  geom_smooth(method = "lm", aes(color = "Best fit"), se = FALSE) +
  labs(
    x = "Observed",
    y = "Predicted mean",
    color = "Linear model",
    title = paste0("Predicted - Observed: ", ylab),
    subtitle = paste0("R sq.: ", rsq)
  ) +
  my_theme()

dest <- paste0(ylab, "_predictedObserved.jpeg")
save_gg(dest, g, "Plots/GAMs")

plot_slopes(
  fit,
  variables = "I_property_area_x_total_take_density",
  condition = "property_area"
) +
  geom_hline(yintercept = 0, linetype = 3)


avg_slps <- avg_slopes(fit)

slps <- slopes(fit, variables = c("property_area"))
head(slps)

# counterfactual predictions

data$y |> summary()

# Calculating the Sturges bins
breaks <- pretty(range(exp(data$y)), n = nclass.Sturges(exp(data$y)), min.n = 1)
y <- exp(data$y)

hist_y <- density(y, from = min(y), to = max(y)) %$%
  data.frame(x = x, y = y) %>%
  mutate(
    id = if_else(x < 15, "High", "Medium"),
    id = if_else(x > 75, "Low", id),
    id = factor(id, levels = c("Low", "Medium", "High"))
  )

breaks <- seq(0, 150, by = 15)

gg <- hist_y |>
  ggplot(aes(x = x, ymin = 0, ymax = y, fill = id)) +
  geom_ribbon() +
  geom_line(aes(y = y)) +
  labs(
    x = "Percent error",
    y = "Density",
    title = "Percent error classifications",
    fill = "Confidence"
  ) +
  scale_fill_manual(
    values = c("Low" = "#fc8d59", "Medium" = "#ffffbf", "High" = "#91bfdb")
  ) +
  coord_cartesian(xlim = c(0, 150)) +
  scale_x_continuous(breaks = breaks) +
  my_theme(16)
save_gg("mpe_confidence", gg, "Plots")

data |>

  ggplot(densities, aes(x = x, y = y, group = id)) +
  geom_density(stat = 'identity') +
  geom_density(
    aes(fill = group),
    . %>%
      filter(
        (group == "C" & between(x, -1.5, 2.0)) |
          (group == "P" & between(x, 0.5, 2.8))
      ),
    stat = 'identity',
    alpha = 0.75
  )

# coef <- coefficients(fit)
# slopes <- slopes(fit)
# slopes_avg <- avg_slopes(fit)
# slopes_avg

# overdisp_fun <- function(model) {
#   rdf <- df.residual(model)
#   rp <- residuals(model,type="pearson")
#   Pearson.chisq <- sum(rp^2)
#   prat <- Pearson.chisq/rdf
#   pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
#   c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
# }
#
# overdisp_fun(fit)
#
# quasi_table <- function(model,ctab=coef(summary(model)),
#                         phi=overdisp_fun(model)["ratio"]) {
#   qctab <- within(as.data.frame(ctab),
#                   {   `Std. Error` <- `Std. Error`*sqrt(phi)
#                   `z value` <- Estimate/`Std. Error`
#                   `Pr(>|z|)` <- 2*pnorm(abs(`z value`), lower.tail=FALSE)
#                   })
#   return(qctab)
# }
#
# quasi_table(fit) |>
#   mutate(sig = `Pr(>|z|)` < 0.05)

# Compute adjusted predictions for each observed combination of regressor in the dataset
# used to fit the model
