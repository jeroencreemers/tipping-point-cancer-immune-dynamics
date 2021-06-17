### Figure 1 ###
# Libraries
library(tidyverse)

# Function to call CPP model from R
source("call_model_function.R")

df_clearance <- call.model(R = 1, xi = 0.005)
df_control <- call.model(R = 1, xi = 0.0005)
df_outgrowth <- call.model(R = 1, xi = 0.00025)

# Bind dataframes + create column with id
df <- bind_rows( list( Clearance = df_clearance,
                       Outgrowth = df_outgrowth,
                       Control =df_control ),
                 .id = "id" ) %>%
  
  # Wide to long format
  gather(`tumor_cells`,`immune_cells`, `specific_cells`, `naive_cells`, key = "cell_type", value = "cells")

# Save dataframe
save(list = "df", file = "fig_1_data.RData")
