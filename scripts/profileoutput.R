file_path <- "profiler/input/Example1.inp"
lines <- readLines(file_path)

# Skip header lines (e.g., assume first 15 lines are headers)
data_lines <- lines[16:length(lines)]
data_lines <- trimws(data_lines)
# Split each line into individual numbers (assuming whitespace-delimited)
data_list <- strsplit(data_lines, "\\s+")

# Convert the list into a matrix, then a data frame
# Use `do.call()` to rbind the list of split lines
data_matrix <- do.call(rbind, data_list)

# Convert matrix to data frame and specify column names
data_frame <- as.data.frame(data_matrix, stringsAsFactors = FALSE)

# Assign column names (adjust these names based on your data)
colnames(data_frame) <- c("Depth", "Porosity", "Biodiffusivity", "Irrigation", "Concentration")

# Convert columns to numeric (if they are numeric data)
data_frame[] <- lapply(data_frame, as.numeric)

# View the resulting data frame
print(data_frame)