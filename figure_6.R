library(survival)
library(MASS)
library(dplyr)
library(ggplot2)
library(survminer)
library(patchwork)

load("figure_6_data.RData")

x_axis_breaks <- c(200, 400, 800, 1600)

layers_for_ab <- list(geom_point(size = 0.75, col = 1+data$long_survival),
                        scale_x_continuous(name = "LDH (U/L)", 
                                           trans = "log", 
                                           limits = c(200, 2000), 
                                           breaks = x_axis_breaks),
                        scale_y_continuous(name = "I/P ratio", trans = "log", 
                                           breaks = c(0.01, 0.1, 1, 10, 100), 
                                           labels = c(0.01, 0.1, 1, 10, 100)))

### Plots A-C ###
plot_6a <- ggplot(data = data, aes(x = LDH, y = Ratio)) +
  layers_for_ab +
  geom_vline(xintercept = ldh.cutoff, col = "grey") +
  annotate("text", x = 1250, y = 15, label = accs( table( m1.ldh, data$long_survival ) ))

plot_6b <- ggplot(data = data, aes(x = LDH, y = Ratio)) +
  layers_for_ab +
  geom_hline(yintercept = ratio.cutoff, col = "grey") +
  annotate("text", x = 1250, y = 15, label = accs( table( m1.ratio, data$long_survival ) ) )

plot_6c <- ggplot(data = data, aes(x = log(LDH), y = log(Ratio))) +
  geom_point(size = 0.75, col = 1+data$long_survival) +
  scale_x_continuous(name = "LDH (U/L)", limits = log(c(200, 2000)), breaks = log(x_axis_breaks), labels = x_axis_breaks) +
  scale_y_continuous(name = "I/P ratio", breaks = log(c(0.01, 0.1, 1, 10, 100)), labels = c(0.01, 0.1, 1, 10, 100)) +
  geom_abline(intercept = (t(m.combi$scaling) %*% colMeans( m.combi$means ))/m.combi$scaling[1], 
              slope = ((-1*m.combi$scaling[2])/ m.combi$scaling[1]), col = "grey") +
  annotate("text", x = log(1250), y = log(15), label = accs( table( mc.combi, data$long_survival ) ) )



### Plots D-F ###
# Kaplan meier curves
layers_KM <- list(scale_x_continuous(name = "Time (months)", breaks = seq(0, 48, 6), limits = c(0,48)))

# Plot 6d
surv_obj <- Surv( data$survivalMonths, data$survivalStatus )
m <- predict(coxph( surv_obj ~ log(LDH), data ))
fit <- survfit( surv_obj ~ cutl( m, 2 ), data = data)

plot_6d <- ggsurvplot(fit = fit, data = data, size = 0.5, xlim = c(0, 48), 
                      palette = c("red", "black"))$plot +
  annotate("text", x=36, y = 0.85, size = 2.5, label = bics(coxph( surv_obj ~ log(LDH), data ))) +
  layers_KM

# Plot 6e
m <- predict(coxph( surv_obj ~ log(Ratio), data ))
fit <- survfit( surv_obj ~ cutl( m, 2 ), data = data)

plot_6e <- ggsurvplot(fit = fit, data = data, size = 0.5, xlim = c(0, 48), 
                      palette = c("red", "black"))$plot +
  annotate("text", x=36, y = 0.85, size = 2.5, label = bics(coxph( surv_obj ~ log(Ratio), data ))) +
  layers_KM

# Plot 6f
m <- predict(coxph( surv_obj ~ log(Ratio)+log(LDH), data ))
fit <- survfit( surv_obj ~ cutl( m, 2 ), data = data)

plot_6f <- ggsurvplot(fit = fit, data = data, size = 0.5, xlim = c(0, 48), 
                      palette = c("red", "black"))$plot +
  annotate("text", x=36, y = 0.85, size = 2.5, label = bics(coxph( surv_obj ~ log(Ratio)+log(LDH), data ))) +
  layers_KM

### Figure 6 ###
figure_6 <- (plot_6a + plot_6b + plot_6c) / (plot_6d + plot_6e + plot_6f) + 
              plot_annotation(tag_levels = 'A')



