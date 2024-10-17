library(ggplot2)
library(reshape2)
library(plotly)
setwd("~/Desktop/CMIP6")
all_files <- list.files(pattern = '*.dat')
string_part1 <- sapply(all_files, function(x) {return (strsplit(x, '_')[[1]][1])})
string_part2 <- sapply(all_files, function(x) {return (strsplit(x, '_')[[1]][2])})
string_part3 <- sapply(all_files, function(x) {return (strsplit(x, '_')[[1]][3])})
string_part4 <- sapply(all_files, function(x) {return (strsplit(x, '_')[[1]][4])})
string_part5 <- sapply(all_files, function(x) {return (strsplit(x, '_')[[1]][5])})
string_part6 <- sapply(all_files, function(x) {return (strsplit(x, '_')[[1]][6])})
string_part7 <- sapply(all_files, function(x) {return (strsplit(x, '_')[[1]][7])})
all_files_split <- data.frame(v1 = string_part1, v2 = string_part2, v3 = string_part3,
                              v4 = string_part4, v5 = string_part5, v6 = string_part6, v7 = string_part7)

# real_data <- read.table('had4_krig_annual_v2_0_0.txt')
# real_data$V2 <- real_data$V2 + 287.75 - real_data$V2[156]
# real_data$V6 <- real_data$V2 - 273.15
# save(real_data, file = 'real_data.Rda')


all_data <- list()
count <- 0
for (file in all_files)
{
  temp <- read.table(file, header = F)
  count <- count + 1
  if (all_files_split$v5[count] == 'historicalNat') {next}
  all_data[[paste0('_s_', all_files_split$v4[count], '_s_', all_files_split$v5[count], 
                   '_s_', all_files_split$v6[count])]] <- temp
  print(dim(temp))
  
}



max_year <- 0
min_year <- 2000
for (data_set in all_data)
{
  max_year <- max(max_year, data_set$V1)
  min_year <- min(min_year, data_set$V1)
}

cleaned_all_data <- list()
months_average <- as.character(months(as.Date('2017-01-01') + 31 * 0:11))
months_average[13] <- 'Average'
for (i in 1:12)
{
  cleaned_all_data[[months_average[i]]] <- data.frame(year = min_year:max_year)
  count <- 2
  for (model in names(all_data))
  {
    cleaned_all_data[[months_average[i]]] <- merge(cleaned_all_data[[months_average[i]]], 
                                                    all_data[[model]][, c('V1', paste0('V', i+1))], #paste0('V', i+1) means the month
                                                    by.x = 'year', 
                                                    by.y = 'V1', all.x = TRUE)
    names(cleaned_all_data[[months_average[i]]])[count] <- model
    if (any(cleaned_all_data[[months_average[i]]][, count] < 273, na.rm = T))
    {
      cleaned_all_data[[months_average[i]]][which(cleaned_all_data[[months_average[i]]][, count] < 273), count] <- NA
    }
    count <- count + 1
  }
  cleaned_all_data[[months_average[i]]]$Quant025 <- apply(cleaned_all_data[[months_average[i]]][,-1], 1, quantile, probs = 0.025, na.rm = TRUE)
  cleaned_all_data[[months_average[i]]]$Quant050 <- apply(cleaned_all_data[[months_average[i]]][,-1], 1, quantile, probs = 0.05, na.rm = TRUE)
  cleaned_all_data[[months_average[i]]]$Quant500 <- apply(cleaned_all_data[[months_average[i]]][,-1], 1, quantile, probs = 0.5, na.rm = TRUE)
  cleaned_all_data[[months_average[i]]]$Quant950 <- apply(cleaned_all_data[[months_average[i]]][,-1], 1, quantile, probs = 0.95, na.rm = TRUE)
  cleaned_all_data[[months_average[i]]]$Quant975 <- apply(cleaned_all_data[[months_average[i]]][,-1], 1, quantile, probs = 0.975, na.rm = TRUE)
  cat(months_average[i], '\n')
}

cleaned_all_data$Average <- cleaned_all_data[[1]][, 1:1322]# without quantile, Jan
for (i in 2:12)
{
  cleaned_all_data$Average <- cleaned_all_data$Average + cleaned_all_data[[i]][, 1:1322]
}
cleaned_all_data$Average <- cleaned_all_data$Average / 12
cleaned_all_data$Average$Quant025 <- apply(cleaned_all_data$Average[,-1], 1, quantile, probs = 0.025, na.rm = TRUE)
cleaned_all_data$Average$Quant050 <- apply(cleaned_all_data$Average[,-1], 1, quantile, probs = 0.05, na.rm = TRUE)
cleaned_all_data$Average$Quant500 <- apply(cleaned_all_data$Average[,-1], 1, quantile, probs = 0.5, na.rm = TRUE)
cleaned_all_data$Average$Quant950 <- apply(cleaned_all_data$Average[,-1], 1, quantile, probs = 0.95, na.rm = TRUE)
cleaned_all_data$Average$Quant975 <- apply(cleaned_all_data$Average[,-1], 1, quantile, probs = 0.975, na.rm = TRUE)


### calculate Average_Diff
a <- cleaned_all_data$Average[, 2:1322]
column_names <- colnames(a)
selected_model <- sub("^_s_(.*?)_s_.*$", "\\1", column_names)
unique_elements <- unique(selected_model)
unique_indices <- match(unique_elements, selected_model)
unique_indices <- unique_indices
mean_results <- list()
for (i in seq_along(unique_indices)) {
  start_index <- unique_indices[i]
  if (i < length(unique_indices)) {
    end_index <- unique_indices[i + 1] - 1
  } else {
    end_index <- ncol(a)  # Extend to the last column if this is the last index
  }
  
  # Compute the mean across the specified columns
  column_means <- rowMeans(a[, start_index:end_index, drop = FALSE])
  
  # Store the result
  mean_results[[i]] <- column_means - 273.15
}
result_dataframe <- as.data.frame(do.call(cbind, mean_results))
result_dataframe <- cbind(cleaned_all_data$Average[, 1], result_dataframe)
names(result_dataframe) <- c("year", unique_elements) # totally 54 model, with one special

cleaned_all_data$Average_Diff <- result_dataframe[1:156,] # year to 2005
cleaned_all_data$Average_Diff$Quant025 <- apply(cleaned_all_data$Average_Diff[,-1], 1, quantile, probs = 0.025, na.rm = TRUE)
cleaned_all_data$Average_Diff$Quant050 <- apply(cleaned_all_data$Average_Diff[,-1], 1, quantile, probs = 0.05, na.rm = TRUE)
cleaned_all_data$Average_Diff$Quant500 <- apply(cleaned_all_data$Average_Diff[,-1], 1, quantile, probs = 0.5, na.rm = TRUE)
cleaned_all_data$Average_Diff$Quant950 <- apply(cleaned_all_data$Average_Diff[,-1], 1, quantile, probs = 0.95, na.rm = TRUE)
cleaned_all_data$Average_Diff$Quant975 <- apply(cleaned_all_data$Average_Diff[,-1], 1, quantile, probs = 0.975, na.rm = TRUE)
cleaned_all_data$Average_Diff$mean <- rowMeans(cleaned_all_data$Average_Diff[, -1], na.rm = TRUE)

cleaned_all_data_2024 = cleaned_all_data

path <- "~/Documents/UW Courses/Research/CO2_Data/cleaned_all_data_2024.Rda"
save(cleaned_all_data_2024, file = path)




