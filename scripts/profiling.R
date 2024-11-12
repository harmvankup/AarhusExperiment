library(tidyverse)
library(ggplot2)
library(readxl)

# Define data inputs (replace these paths with your actual file paths)
con_vs_depth_file <- read_xlsx("input/ConVsDepth.xlsx", col_names = F)

porosity_vs_depth_file <- read.csv("input/PorVsDepth.csv")


# Preset values
Datatype <- 2  # Set data type (0, 1, 2, 3 corresponding to choices in Shiny app)
diffusivity_change_with_depth <- FALSE
compaction_considered <- FALSE
FF <- 2.2
Por <- 0.9
Diff <- 0.018
Temp <- 20
sedimentation_rate <- 0
porosity_reference <- 0
flow <- 0
significance_level <- 0.05
zone_number <- 10
mgn <- 10

# Load data
con_vs_depth <-con_vs_depth_file
depth_data <- as.matrix(con_vs_depth)[, 1]
con_data <- as.matrix(con_vs_depth)[, 2]

# Load formation factor and porosity data if available
formation_factor_vs_depth <-  NULL
porosity_vs_depth <- if (Datatype %in% c(1, 2)) as.matrix(porosity_vs_depth_file) else NULL
# Process data based on input type
ff_depth <- porosity_vs_depth[, 1]
ff_data <- 10^(.1916) * porosity_vs_depth[, 2]^(-1.8812)

por_depth <-  porosity_vs_depth[, 1]
por_data <-  porosity_vs_depth[, 2]

# Define functions for calculations
pad_function <- function(data) {
  early_pad <- rev(data[1:3])
  early_pad <- data[1] + (data[1] - early_pad)
  
  late_pad <- rev(data[(length(data) - 2):length(data)])
  late_pad <- data[length(data)] + (data[length(data)] - late_pad)
  
  padded <- c(early_pad[1:2], data, late_pad[2:3])
  smoothed <- padded %>% 
    as_tibble() %>% 
    rename_at(1, ~"concentration") %>% 
    transmute(smoothed = lag(concentration, n = 2) * 0.06 +
                0.24 * lag(concentration, n = 1) +
                0.4 * concentration +
                0.24 * lead(concentration, 1) +
                0.06 * lead(concentration, 2)) %>%
    filter(!is.na(smoothed)) %>% 
    as.matrix()
  return(smoothed)
}


# Smooth data
smoothed_con_data <- pad_function(con_data)
smoothed_ff_data <-  pad_function(ff_data)
smoothed_por_data <-  pad_function(por_data)

# Interpolation and layer calculations
layer_fn <- function(data) {
  layer <- mean(diff(data)) / 2
  while (layer >= 1*mgn) {
    layer <- layer / 2
  }
  return(layer)
}

layer <- layer_fn(depth_data)
interp_depth <- seq(min(depth_data), max(depth_data), by = layer)

interp_fn <- function(depth, data) {
  approx(depth, data, xout = interp_depth, method = "linear", rule = 2)$y
}

interp_con_data <- interp_fn(depth_data, smoothed_con_data)
interp_ff_data <- interp_fn(ff_depth, smoothed_ff_data)
interp_por_data <- interp_fn(por_depth, smoothed_por_data)

# Diffusivity calculations
viscosity_0 <- if (!diffusivity_change_with_depth) {
  1
} else if (0 <= Temp && Temp < 20) {
  10^((1301 / (998.333 + 8.1855 * (Temp - 20) + 0.00585 * (Temp - 20)^2)) - 3.30233)
} else if (20 < Temp && Temp <= 100) {
  1.002 * 10^((1.3272 * (20 - Temp) - 0.001053 * (Temp - 20)^2) / (Temp + 105))
} else if (Temp == 20) {
  1.002
} else {
  stop("Temperature exceeds calculation range (0-100 Celsius)")
}

viscosity <- if (!diffusivity_change_with_depth) {
  rep(1, length(depth_data))
} else {
  approx(depth_data, rep(Temp, length(depth_data)), xout = depth_data, method = "linear", rule = 2)$y %>%
    sapply(function(temp) {
      if (0 <= temp && temp < 20) {
        10^((1301 / (998.333 + 8.1855 * (temp - 20) + 0.00585 * (temp - 20)^2)) - 3.30233)
      } else if (20 < temp && temp <= 100) {
        1.002 * 10^((1.3272 * (20 - temp) - 0.001053 * (temp - 20)^2) / (temp + 105))
      } else if (temp == 20) {
        1.002
      } else {
        stop("Temperature exceeds calculation range (0-100 Celsius)")
      }
    })
}

interp_diffusivity <- if (!diffusivity_change_with_depth) {
  rep(Diff, length(interp_depth))
} else {
  ((Temp + 273.5) / viscosity) * viscosity_0 * Diff / (Temp + 273.5)
}

interpburial <-  if (Datatype == 0 | compaction_considered == F){
  sedimentation_rate
} else {
  porosity_reference * (1 - smoothed_por_data[1]) * sedimentation_rate / (1 - porosity_reference)
}

interpw <- 
  flow* smoothed_por_data[1]



# Burial velocity calculations
burial_velocity <- if (Datatype == 0 || !compaction_considered) {
  rep(sedimentation_rate, length(interp_depth))
} else {
  porosity_reference * (1 - smoothed_por_data[1]) * sedimentation_rate / (1 - porosity_reference)
}

# Reaction rates (simplified example)
reaction_rates <- -interp_con_data * burial_velocity

# Define tridiagonal matrix components
AA <- rep(0, length(interp_depth))
BB <- rep(0, length(interp_depth))
CC <- rep(0, length(interp_depth))

for (i in 2:(length(interp_depth) - 1)) {
  AA[i] <- 2 * interp_diffusivity[i + 1] / interp_ff_data[i + 1] *
    (2 * interp_diffusivity[i] / interp_ff_data[i] - (burial_velocity[i] + flow) * layer) /
    (2 * interp_diffusivity[i + 1] / interp_ff_data[i + 1] * layer +
       2 * interp_diffusivity[i] / interp_ff_data[i] * layer) / layer
  
  BB[i] <- -2 * interp_diffusivity[i] / interp_ff_data[i] *
    (2 * interp_diffusivity[i + 1] / interp_ff_data[i + 1] + (burial_velocity[i + 1] + flow) * layer) /
    (2 * interp_diffusivity[i + 1] / interp_ff_data[i + 1] * layer +
       2 * interp_diffusivity[i] / interp_ff_data[i] * layer) / layer
  
  CC[i] <- 2 * interp_diffusivity[i - 1] / interp_ff_data[i - 1] *
    (2 * interp_diffusivity[i] / interp_ff_data[i] + (burial_velocity[i] + flow) * layer) /
    (2 * interp_diffusivity[i] / interp_ff_data[i] * layer +
       2 * interp_diffusivity[i - 1] / interp_ff_data[i - 1] * layer) / layer
}

# Construct tridiagonal matrix A
A_matrix <- matrix(0, nrow = length(interp_depth) - 1, ncol = length(interp_depth) - 1)
diag(A_matrix) <- BB[2:(length(BB) )]
A_matrix[cbind(2:(length(interp_depth) - 1), 1:(length(interp_depth) - 2))] <- AA[3:(length(AA) - 1)]
A_matrix[cbind(1:(length(interp_depth) - 2), 2:(length(interp_depth) - 1))] <- CC[3:(length(CC) - 1)]

# Helper vector f
f_vec <- rep(0, length(interp_depth) - 1)

# SSE Function
SSErev <- function(dbest, f, interpdepth, interpCondata, AA, BB, CC, index, depthdata, smoothedCondata, Condata) {
  n <- length(interpdepth)
  newd  <- numeric(n)
  newd[1:index[1]] <- dbest[1]
  if (length(dbest) > 1) {
    for (j in 2:length(dbest)) {
      start <- index[j - 1] + 1
      end <- index[j]
      newd[start:end] <- dbest[j]
    }
  }
  
  e <-  rep(0, n - 1)
  e[2] <- AA[2] / BB[2]
  
  for (i in 3:(n - 2)) {
    e[i] <- AA[i] / (BB[i] - e[i - 1] * CC[i])
  }
  
  newd <- as.numeric(newd)
  y <- numeric(n)
  y[1] <- interpCondata[1]
  y[n] <- interpCondata[n]
  
  lf <-  length(f)
  f[2] <- (newd[2] - CC[2] * interpCondata[1]) / BB[2]
  
  for (m in 3:(lf - 1)) {
    f[m] <- (newd[m] - f[m - 1] * CC[m]) / (BB[m] - e[m - 1] * CC[m])
  }
  
  f[lf] <- (newd[lf] - AA[lf] * interpCondata[length(interpCondata)] - f[lf - 1] * CC[lf]) /
    (BB[lf] - e[lf - 1] * CC[lf])
  
  y[lf] <- f[lf]
  
  for (n in (lf - 1):2) {
    y[n] <- f[n] - e[n] * y[n + 1]
  }
  
  # Sum of squared errors
  l <- length(depthdata)
  
  confit <- approx(x = interpdepth, y = y, xout = depthdata, method = "linear", rule = 2, ties = "ordered")$y
  
  sumfit <- sum(confit)
  summeasure <- sum( Condata)
  
  meanfit <- sumfit / l
  meanmeasure <- summeasure / l
  
  err <- numeric(l)
  
  
  for (o in 1:l) {
    err[o] <- (confit[o] - smoothedCondata[o])^2
  }
  SSE1 <- sum(err)
  
  Sxy <- 0
  Sx <- 0
  Sy <- 0
  
  for (p in 1:length(depthdata)) {
    Sxy <- Sxy + (confit[p] - meanfit) * (Condata[p] - meanmeasure)
    Sx <- Sx + (confit[p] - meanfit)^2
    Sy <- Sy + (Condata[p] - meanmeasure)^2
  }
  
  if (Sx != 0 & Sy != 0) {
    R_squared <- Sxy^2 / Sx / Sy
  } else {
    R_squared <- 0
  }
  

  # Assign computed values back to the environment
  assign("y", y, pos = 1)
  assign("newd", newd, pos = 1)
  assign("R_squared", R_squared, pos = 1)
  assign("f",f,pos =1)
  assign("confit",confit, pos = 1)
  assign("e",e, pos = 1)
  return(SSE1)
}

# Optimization loop for zone optimization
Maximum_Number_of_Zones <-  floor(length(depth_data / zone_number))
max_zones <- 5
best_sse <- Inf
best_num_zones <- 1

tol <- 1e-6 * mean(con_data) 

match <- matrix(0, nrow = max_zones, ncol = max_zones)

for (num_zones in 1:max_zones) {
  
  # Split data into zones
  n <- num_zones
  D <- rep(0, num_zones)
  dbest <- D[1:num_zones]
  zone_length <- floor(length(interp_depth) / num_zones)
  zonebreak <- rep(0, nzones)
  for (j in 1:nzones) {
    zonebreak[j] <- depthdata()[j * zonelength]
  }
  zonebreak[num_zones] <- depth_data[length(depth_data)]
  zone_indices <- split(1:length(interp_depth), ceiling(seq_along(1:length(interp_depth)) / zone_length))
  index <- match[1:n, num_zones]
  # Optimization using Nelder-Mead method for each zone
  for (zone in zone_indices) {
    optim_result <- optim(
      par = dbest,
      fn = SSErev,
      f = f_vec,
      index = index,
      interpdepth = interp_depth,
      interpCondata = interp_con_data,
      depthdata = depth_data,
      Condata = con_data,
      smoothedCondata = smoothed_con_data,
      AA = AA,
      BB = BB,
      CC = CC,
      method = "Nelder-Mead"
    )
    
    # Calculate SSE for the current optimization
    current_sse <- optim_result$value
    if (current_sse < best_sse) {
      best_sse <- current_sse
      best_num_zones <- num_zones
    }
  }
}

# Extract optimized values from the best zone configuration
optimized_y <- numeric(length(interp_depth))
for (zone in zone_indices[1:best_num_zones]) {
  optim_result <- optim(
    par = rep(1, length(zone) - 1),
    fn = define_sse,
    f = f_vec[zone],
    interpdepth = interp_depth[zone],
    interpCondata = interp_con_data[zone],
    AA = AA[zone],
    BB = BB[zone],
    CC = CC[zone],
    method = "Nelder-Mead"
  )
  optimized_y[zone] <- optim_result$par
}

# Plot with optimized results
ggplot() +
  geom_line(aes(x = con_data, y = depth_data), color = "blue") +
  geom_point(aes(x = con_data, y = depth_data), color = "red", shape = 1) +
  geom_line(aes(x = reaction_rates, y = interp_depth), color = "green") +
  geom_line(aes(x = optimized_y, y = interp_depth), color = "purple", linetype = "dashed") +
  scale_y_reverse() +
  labs(x = "Concentration / Reaction Rate", y = "Depth", title = "Concentration, Reaction Rate, and Optimized Results vs Depth") +
  theme_minimal()
