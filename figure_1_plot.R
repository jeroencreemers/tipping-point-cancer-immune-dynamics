### Figure 1 - Plot ###
# Libraries
library(ggplot2)
library(scales)

# Load data
load(file = "fig_1_data.RData")

# Plot 
ggplot(df, aes(x=time, y = cells, col = cell_type, linetype = cell_type )) + 
  geom_hline(yintercept = 10^12, linetype = 2, color = "grey", size = 0.2) + # time of diagnosis
  geom_hline(yintercept = 50*10^8, linetype = 2, color = "grey", size = 0.2) + # time of death
  geom_line() +
  facet_grid(~id, margins = c(0,0,0,0)) +
  scale_color_manual(values=c(naive_cells="#019939", specific_cells = "#3c63ad", immune_cells = "#3c63ad", tumor_cells = "#cc6615")) +
  scale_linetype_manual(values=c(naive_cells="dotted", specific_cells = "dotted", immune_cells = "solid", tumor_cells = "solid")) +
  scale_y_continuous(name = "Cells",
                     expand = c(0,0),
                     trans = "log10", 
                     breaks = c(10^0, 10^3, 10^6, 10^9, 10^12), 
                     labels = trans_format("log10", math_format(10^.x)), 
                     limits =c(0.1, 10^13)) +
  scale_x_continuous(name = "Time (months)", expand = c(0,0),  
                     limits = c(0, (5*365)), 
                     breaks = c(0, 365, (2*365), (3*365), (4*365), (5*365)), 
                     labels = c(0, 12, 24, 36, 48, 60)) +
  theme(panel.spacing = unit(1, "lines"), strip.text.x = element_text(size = 9), 
        strip.background = element_blank()) 