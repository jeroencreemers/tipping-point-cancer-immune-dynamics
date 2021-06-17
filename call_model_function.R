# Function to call C++ model from within R 
call.model <- function(R=5, xi=0.005, stochastic_killing=0, stochastic_growth=0,
                       raise_killing=1, treatment_duration=730, drift_r=0, 
                       drift_xi=0, seed = 0){
  read.table(text=system2("../model/./tumormodel",c("--R",R,
                                                    "--xi",xi,
                                                    "--stochastic-killing", stochastic_killing, 
                                                    "--stochastic-growth", stochastic_growth, 
                                                    "--raise-killing", raise_killing, 
                                                    "--treatment-duration", treatment_duration, 
                                                    "--drift-R", drift_r, 
                                                    "--drift-xi", drift_xi, 
                                                    "--seed", seed), 
                          stdout=T)) %>% 
    rename(time=V1, 
           R = V2, 
           xi = V3, 
           tumor_cells = V4, 
           immune_cells = V5, 
           specific_cells = V6, 
           naive_cells = V7) %>%
    
    mutate(status = ifelse(max(tumor_cells > 10^12), yes = 1, no = 0)) %>%
    
    filter(tumor_cells < 10^12) %>%
    
    # Find the first value of time at which tumor cells > 65*10^8 
    mutate(time_diagnosis = time[tumor_cells > 65*10^8][1], 
           OS = max(time) - time_diagnosis)
}