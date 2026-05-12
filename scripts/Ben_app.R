############################################################
## JTB MASTER PIPELINE: AUTONOMOUS SURVEILLANCE OPTIMIZATION
## Case Study: Cassava Brown Streak Disease (CBSD) in Benin
############################################################

############################################################
## JTB PIPELINE: RESOLVED GEOMETRY & AUTONOMOUS SAMPLING
############################################################

library(sf)
library(tidyverse)
library(geodata)
library(viridis)
library(patchwork)

# --- STAGE 1: GEOMETRY RESOLUTION ---
# Fetching and cleaning the Benin map to prevent the "Edge Crosses Edge" error
benin_sf <- gadm(country = "BEN", level = 1, path = tempdir()) %>% 
  st_as_sf() %>% 
  st_make_valid() # Critical fix for S2 engine

# --- STAGE 2: THE ENGS AUTO-SOLVER ---
# n = farms, m = plants. Automatically determined by Budget B=4000
solve_design <- function(theta, phi) {
  expand.grid(n = 5:1000, m = 10:30) %>%
    mutate(Cost = n * 150 + n * m * 5) %>% # B = n*Ct + n*m*Cp
    filter(Cost <= 40000) %>%
    mutate(ENGS = (n * (1 - (1 - phi * theta)^m) * 2000) - Cost) %>%
    filter(ENGS == max(ENGS)) %>% slice(1)
}

# --- STAGE 3: MAPPING THE GEOMETRY OF SAMPLING ---
run_benin_analysis <- function(phase) {
  # Biological state from DDE logic
  state <- if(phase == "Onset") list(t=0.1, p=0.15) else list(t=0.4, p=0.6)
  
  # 1. Get optimal n* and m*
  opt <- solve_design(state$t, state$p)
  
  # 2. Create sampling grid
  # We use st_make_valid again here to ensure the grid-map intersection works
  grid <- st_make_grid(benin_sf, n = 1000) %>% 
    st_as_sf() %>% 
    st_make_valid() %>% 
    st_intersection(st_union(benin_sf))
  
  # 3. Calculate Information Frontier (W_j)
  coords <- st_coordinates(st_centroid(grid))
  # Onset focuses on the Eastern Border (Nigeria)
  mu <- if(phase == "Onset") 0.15 * exp(-0.05 * (max(coords[,1]) - coords[,1])) else 0.4
  grid$W_j <- mu * (1 - mu)
  
  # 4. Deploy n* farms
  set.seed(123)
  points <- grid[sample(1:nrow(grid), size = opt$n, prob = grid$W_j), ]
  
  # 5. Visual Result
  ggplot() +
    geom_sf(data = grid, aes(fill = W_j), color = NA) +
    geom_sf(data = points, color = "white", size = 2, shape = 21, fill = "black", stroke = 1) +
    scale_fill_viridis_c(option = "mako", name = "Wj") +
    labs(title = paste(phase, "Phase"),
         subtitle = paste("n* =", opt$n, "farms | m* =", opt$m, "plants/farm")) +
    theme_void()
}

# EXECUTE
(run_benin_analysis("Onset") | run_benin_analysis("Established"))