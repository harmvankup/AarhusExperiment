# Define header lines exactly as seen in the example file
header <- c(
  "Interpretation of O2-profile in paper (Fig. 3)",
  "-0.02     Depth at top of calculation domain",
  " 0.27     Depth at bottom of calculation domain",
  " 7   Max number of equally spaced zones in interpretation (1 to 12)",
  " 3        Type of boundary conditions (1:t=C b=C, 2:t=C t=F, 3:b=C b=F 4:t=C b=F 5:t=F b=C)",
  " 0.22573  First boundary condition",
  " 0.0      Second boundary condition",
  "11.7E-06  Diffusivity in water (D)",
  " 2        Expression for sediment diffusivity (Ds) (1: Ds=FI*D, 2: Ds=FI**2*D, 3: Ds=D/(1+3*(1-FI))",
  "321.4123  Concentration in water column (C0)",
  "-1.0E+20  Minimum for production rate",
  " 0.0      Maximum for production rate",
  " 0.001    Maximum deviation (in %) when accepting a calculated minimum",
  " 0.01     Level of significance in the F statistic",
  "    X    FI   DB  ALFA           C"
)

# Example data frame for the data section
depth_concentration <- data_frame

# Format the data frame to match the fixed-width structure
formatted_data <- apply(depth_concentration, 1, function(row) {
  sprintf("%7.2f %5.0f %4.0f %5.0f %13.4f", row[1], row[2], row[3], row[4], row[5])
})

# Open a connection to the .inp file
file_connection <- file("profiler/input/test.inp", "w")

# Write the header
writeLines(header, file_connection)

# Write the formatted data
writeLines(formatted_data, file_connection)

# Close the connection
close(file_connection)

cat("File 'test.inp' has been created.\n")
