library(tidyverse)
library(ggforce)
library(ggprism)
library(stringr)
library(ggplot2)
library(grid)
library(gridExtra)
library(qqman)
library(ggrepel)
library(optparse)

# List the list of options that the script can accept

option_list <- list(
  make_option(c("--input-directory"), type = "character", default = ".",
              help = "Input directory", metavar = "character"),
  make_option(c("--output-directory"), type = "character", default = ".",
              help = "Output directory", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

print("Options:")
print(opt)

input_directory <- opt$`input-directory`
output_directory <- opt$`output-directory`

print("Input directory:")
print(input_directory)

print("Output directory:")
print(output_directory)

if (!file.exists(input_directory)) {
  stop("Input directory does not exist.")
}

if (!file.exists(output_directory)) {
  print("Output directory does not exist. Creating it.")
  dir.create(output_directory, showWarnings = FALSE, recursive = TRUE)
}

if (grepl("Beijing", input_directory, ignore.case = TRUE)) {
  region <- "Beijing"
} else if (grepl("Guangzhou", input_directory, ignore.case = TRUE)) {
  region <- "Guangzhou"
} else {
  region <- "Unknown"
}

print("Region:")
print(region)

input_files <- list.files(input_directory, pattern = "adjusted$", full.names = TRUE) # nolint
output_file <- ifelse(region == "Beijing",
                      file.path(output_directory, "bac_age_BJ.pdf"),
                      ifelse(region == "Guangzhou",
                             file.path(output_directory, "bac_age_GZ.pdf"),
                             file.path(output_directory, "bac_age.pdf")))

output_file_png <- ifelse(region == "Beijing",
                          file.path(output_directory, "bac_age_BJ.png"),
                          ifelse(region == "Guangzhou",
                                 file.path(output_directory, "bac_age_GZ.png"),
                                 file.path(output_directory, "bac_age.png")))

print("Input files:")
print(input_files)

print("Output file:")
print(output_file)

known_columns <- c("CHR", "SNP", "UNADJ", "GC", "BONF", "HOLM", "SIDAK_SS", "SIDAK_SD", "FDR_BH", "FDR_BY", "Bacteria") # nolint

df_list <- lapply(input_files, function(file) {
  print(paste("Processing file", file))
  df <- read.table(file, header = TRUE, row.names = NULL, sep = ",")
  df <- df[, known_columns[known_columns %in% colnames(df)]]
  df$SNP <- as.character(df$SNP)  # Ensure SNP column is character
  print(paste(file, "has", nrow(df), "rows and", ncol(df), "columns"))
  print(paste("The", file, "has been processed."))
  print(" ")
  return(df)
})

print("Head of the processed data:")
print(head(df_list[[1]]))

print("Processing data...")

df_all <- bind_rows(df_list)

df_all$BP <- sapply(df_all$SNP, function(x) {
  parts <- strsplit(x, ":")
  return(parts[[1]][2])  # 返回分割后的第二部分
})

print("Filtering data...")

# Filter out SNPs with -log10(BONF) > 10, as BONF < 1e-10
df_all <- df_all %>%
  filter(df_all$BONF < 1e-10)

print("Sorting data...")

df_all$CHR[!df_all$CHR %in% c(1:22, 23, 26)] <- 30

df_all$CHR<-as.numeric(df_all$CHR)
df_all$BP<-as.numeric(df_all$BP)

sorted_df_all <- df_all %>%
  arrange(CHR, BP)
sorted_df <- sorted_df_all[, c("CHR", "BP", "SNP", "BONF", "Bacteria")]

name <- unique(sorted_df_all$CHR)

# 将CHR中的23替换为X, 26替换为Y, 其他替换为Others
# old_elements <- c("23", "26", "30") # nolint
# new_elements <- c("X", "Y", "Others") # nolint

#定义df
df <- sorted_df_all

# df$CHR <- ifelse(df$CHR %in% old_elements, new_elements, df$CHR) # nolint

print("Calculating rank...")

df$Bacteria <- as.character(df$Bacteria)

df <- df %>%
  group_by(Bacteria) %>%
  mutate(rank_p = rank(BONF, ties.method = "first")) %>% # 对P进行排名
  mutate(Bacteria_new = if_else(rank_p == 1, Bacteria, NA_character_)) %>% # 保留P最小的行的Bacteria值，其他变成NA # nolint
  ungroup() %>%
  select(-rank_p)

print("Drawing plot...")

# Manhattan Plot的绘制

chromosome_colors <- c(
  rep(c('#1F77B4','#FF7F0C','#2BA099','#D62628','#9467BD','#8C564B','#7F7F7F','#E477C2','#BDBD21', '#17BECF'),3)) # nolint

significant_snps <- df %>%
  filter(!is.na(Bacteria_new)) %>%
  select(CHR, SNP ,BP, Bacteria_new)

snpsOfInterest<-significant_snps$SNP

p <- manhattan(df,
  chr = "CHR",
  bp = "BP",
  p = "BONF",
  snp = "SNP",
  col = chromosome_colors,
  highlight = snpsOfInterest
) 

myplot<-recordPlot()
replayPlot(myplot)

#保存图像
save_manhattan_plot<-function(output_file_pdf,output_file_png){
  pdf(file = output_file_pdf,width = 10,height = 6)
  replayPlot(myplot)
  dev.off()
  png(filename = output_file_png,width = 10*300,height = 6*300,res = 300)
  replayPlot(myplot)
  dev.off()
  cat("Plots saved to:\n",output_file_pdf,"\n",output_file_png,"\n")
}

save_manhattan_plot(output_file,output_file_png)


print("Done.")