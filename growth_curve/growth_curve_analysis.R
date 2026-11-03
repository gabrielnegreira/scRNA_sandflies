#libraries####
library(tidyverse)
library(readODS)
library(ggalign)
library(ggridges)
library(ggrepel)
library(scales)
library(ggh4x)

#source common variables
source("color_pallets.R")

#parameters####
inputs_dir <- "../inputs/" #relative to repo root
outputs_dir <- "../outputs/growth_curves"
statistics_dir <- paste0(inputs_dir, "/growth_curves/statistics/")
wells_metadata <-  paste0(inputs_dir, "/growth_curves/wells_metadata.ods")
cols_to_plot <- c("Has DNA % Parent", 
                  "Is alive % Parent",
                  "Is alive % Total",
                  "Singlets % Parent", 
                  "relative_cell_size", 
                  "relative_DNA_size", 
                  "multiplate_rate",
                  "Singlets Median B610-A", 
                  "Singlets GeoMean FSC-A",
                  "Singlets GeoMean B610-A", 
                  "Is alive SD Time",
                  "original_concentration(par/ml)")

#load files####
##load wells metadata
wells_metadata <- read_ods(wells_metadata, sheet = 1)

## load statistics files
file_names <- list.files(statistics_dir, pattern = "*.csv")
files <- paste0(statistics_dir, file_names)
files <- lapply(files, utils::read.csv, skip = 2, header = FALSE) %>%
  lapply(
    function(x){
      colnames(x) <- x[1,] 
      x <- x[-1,]
      return(x)}
  )
names(files) <- file_names

#bind all files into a single dataframe####
to_plot <- bind_rows(files, .id = "file")

#format the dataframe####
colnames(to_plot) <- gsub(":", "", colnames(to_plot))
to_plot[to_plot == "####"] <- NA
to_plot <- utils::type.convert(to_plot, as.is = TRUE, dec = ",")
##convert percentages to numeric
percent_values <- to_plot[grep("%", to_plot)]
percent_values <- apply(percent_values, 2, function(y){
  y <- gsub("%", "", y)
  #since the table is european-formated, we need to deal with commas and points.
  y <- gsub("\\.", "", y)
  y <- gsub(",", ".", y)
  y <- as.numeric(y)
  y <- y/100
})
to_plot[grep("%", to_plot)] <- percent_values
rm(percent_values)

to_plot <- to_plot %>%
  mutate(
    failed = grepl("_failed", `Tube Name`),
    `Tube Name` = gsub("_failed", "", `Tube Name`),
    `Record Time` = as.POSIXct(`Record Time`, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    day = as.Date(to_plot$`Record Time`, tz = "UTC"),
    day = as.integer(difftime(day, min(day, na.rm = TRUE), units = "days")),
    well = gsub("01-Well-","",`Tube Name`),
    row = str_extract(well, "[A-Z]"),
    row = factor(row, LETTERS[1:8]),
    column = as.numeric(str_extract(well, "\\d+$")),
    column = factor(column, 1:12),
    multiplate_rate = 1 - `Singlets % Parent`
  ) %>%
  group_by(file) %>%
  mutate(relative_cell_size = `Singlets Median FSC-A` / median(`Singlets Median FSC-A`, na.rm = TRUE),
         relative_DNA_size = `Singlets Median B610-A` / median(`Singlets Median B610-A`, na.rm = TRUE))

#append the wells metadata
to_plot <- to_plot %>%
  left_join(wells_metadata, by = "well") %>%
  mutate(subclone = factor(subclone, levels = c("D6", "G8", "E4", "G3", "A9", "G6", "mother clone")))

#plot wells metadata 
plot_list <- list()
for(col in colnames(wells_metadata)){
  if(col == "well"){
    next
  }
  
  plot <- to_plot[,c("row", "column", "file", col)]
  plot_list[[col]] <- plot %>%
    mutate(file = factor(file)) %>%
    filter(as.numeric(file) == 1) %>%
    tidyr::complete(row, column) %>%
    mutate(row = fct_rev(row)) %>%
    ggplot(aes(x = column, y = row, fill = factor(.data[[col]]), label = factor(.data[[col]])))+
    geom_point(size = 20, shape = 21)+
    geom_text()+
    scale_y_discrete(name = NULL, drop = FALSE)+
    scale_x_discrete(name = NULL, position = "top", drop = FALSE)+
    scale_fill_discrete(name = NULL, na.value = "white")+
    ggtitle(col)+
    guides(fill = "none")+
    theme_bw()+
    theme(axis.line = element_blank(), axis.ticks = element_blank())
}
final_plot <- ggalign::align_plots(!!!plot_list)
scale_factor <- 1.5
ggsave(paste0(outputs_dir, "/plate_design.pdf"), final_plot, width = 11.69*scale_factor, height = 8.27*scale_factor, create.dir = TRUE)

#plot the columns specified above
plot_list <- list()
for(col in cols_to_plot){
  plot <- to_plot %>%
    select(row, column, file, all_of(col)) %>%
    mutate(label = .data[[col]])
  
  #format label column
  if(grepl("%", col)){
    plot$label <- paste0(round(plot$label * 100, 1), "%")
  }else{
    plot$label <- round(plot$label, 2)
  }
  
  if(max(plot[[col]], na.rm = TRUE) > 1000){
    plot$label <- format(plot$label, scientific = TRUE, digits = 2)
  }
  
  plot <- plot %>%
    mutate(row = fct_rev(row)) %>%
    group_by(file) %>%
    complete(row, column) %>%
    ggplot(aes(x = column, y = row, fill = .data[[col]], label = label))+
    geom_point(size = 4, shape = 21)+
    geom_text(size = 0.85)+
    scale_y_discrete(name = NULL, drop = FALSE)+
    scale_x_discrete(name = NULL, position = "top", drop = FALSE)+
    scale_fill_viridis_c(name = NULL, option = "turbo", na.value = "white")+
    ggtitle(col)+
    theme_bw()+
    theme(axis.line = element_blank(), axis.ticks = element_blank())
  
  plot_list[[col]] <- plot + facet_grid(cols = vars(file))
  scale_factor <- 1.25
  ncols <- ceiling(sqrt(length(unique(to_plot$file))))
  col <- gsub("%", "percent", col)
  col <- gsub("\\s*\\([^)]*\\)", "", col)
  ggsave(paste0(outputs_dir, "/", col, ".pdf"), plot + facet_wrap(vars(file)), width = 4+(1*ncols*scale_factor), height = 1*ncols*scale_factor, create.dir = TRUE)
}

final_plot <- ggalign::align_plots(!!!plot_list, ncol = 1)
scale_factor <- 2
ggsave(paste0(outputs_dir, "/numeric_variables.pdf"), final_plot, width = 11.69*scale_factor, height = 1*length(cols_to_plot)*scale_factor, create.dir = TRUE)

#estimate average counts per sample
to_plot <- to_plot %>%  
  mutate(original_count = `original_concentration(par/ml)`) %>%
  group_by(sample, day) %>%
  mutate(corrected_count = ifelse(failed, NA, original_count)) %>%
  mutate(corrected_count = ifelse(failed, mean(corrected_count, na.rm = TRUE), corrected_count)) %>%
  mutate(mean_original_count = mean(original_count, na.rm = TRUE)) %>%
  mutate(mean_corrected_count = mean(corrected_count, na.rm = TRUE)) %>%
  ungroup()

#visualize the raw growth curve with the means and failed data points
plot <- to_plot %>%
  filter(day != 5) %>% #day 5 reading failed
  filter(!sample %in% c("blank", "not_stained")) %>%
  filter(! is.na(subclone)) 

plot <- plot %>%
  ggplot(aes(x = as.numeric(day), y = original_count, group = well))+
  geom_line(linewidth = 0.2, color = "lightgrey")+
  geom_line(aes(y = mean_corrected_count), linetype = 'dashed', color = "black", linewidth = 1)+
  geom_text(aes(label = well, color = failed))+
  scale_x_continuous(breaks = unique(plot$day))+
  ggtitle("Raw growth curves", subtitle = "dashed lines: mean between triplicates (failed datapoints ignored); datapoint labels = well in the growth curve plate")+
  facet_wrap(vars(sample))+
  theme_bw()+
  theme(panel.grid.minor.x = element_blank())

ggsave(paste0(outputs_dir, "/raw_growth_curves.pdf"), plot, width = 11.69, height = 9, scale = 2, create.dir = TRUE)

#visualize growth profile of passage 2 parasites per aneuploidy group
plot_data <- to_plot %>%
  filter(subclone != "mother clone") %>%
  filter(day != 5) %>%
  filter(day < max(day, na.rm = TRUE)) %>%
  #filter(day < 7) %>%
  filter(!sample %in% c("blank", "not_stained")) %>%
  filter(! is.na(subclone)) %>%
  select(day, mean_corrected_count, aneu_group, well, subclone, passage, sample) %>%
  distinct() 

label_data <- plot_data %>%
  filter(passage == 2) %>%
  group_by(sample) %>%
  filter(day == 8) %>%
  slice(1) %>%  # In case of ties, take first occurrence
  ungroup()

plot <- plot_data %>%
  filter(passage == 2) %>%
  ggplot(aes(x = day, y = mean_corrected_count, color = aneu_group, group = sample, label = subclone))+
  geom_line()+
  geom_point()+
  geom_text_repel(
    data = label_data, 
    min.segment.length = 0,
    force = 20,
    nudge_x = 2,
    nudge_y = 3,
    color = "black"
  )+
  scale_x_continuous(breaks = unique(plot_data$day))+
  scale_color_manual(values = aneu_colors)+
  facet_grid(rows = vars(passage))+
  labs(y = "live parasites/ml")+
  theme_bw()+
  theme(panel.grid.minor.x = element_blank())

#visualize how each subclone progresses over passages
strip_colors <- plot_data %>%
  filter(passage == 2) %>%
  select(subclone, aneu_group) %>%
  distinct() %>%
  mutate(color = aneu_colors[aneu_group]) %>%
  select(subclone, color) %>%
  deframe()

strip_colors <- strip_colors[levels(plot_data$subclone)]
label_data <- plot_data %>%
  group_by(sample) %>%
  filter(day == 8) %>%
  slice(1) %>%  # In case of ties, take first occurrence
  ungroup()

plot_2 <- plot_data %>%
  ggplot(aes(x = day, y = mean_corrected_count, color = aneu_group, group = sample, label = paste0("P", passage)))+
  geom_line()+
  geom_point()+
  geom_text_repel(
    data = label_data, 
    min.segment.length = 0,
    force = 20,
    nudge_x = 2,
    nudge_y = 3,
    color = "black"
  )+
  facet_wrap2(vars(subclone), nrow = 1, strip = strip_themed(background_x = elem_list_rect(fill = strip_colors))) +
  scale_x_continuous(breaks = unique(plot_data$day))+
  scale_color_manual(values = aneu_colors)+
  labs(y = "live parasites/ml")+
  theme_bw()+
  theme(panel.grid.minor.x = element_blank())

final_plot <- ggalign::align_plots(plot, plot_2, nrow = 1, widths = c(0.2, 0.8))
final_plot
ggsave(paste0(outputs_dir, "_aneuploidy.pdf"), final_plot, width = 11.69, height = 2, scale = 2, create.dir = TRUE)

#plot doubling time
doubling_time <- to_plot %>%
  filter(day <= 2 & day >= 1) %>%
  filter(!sample %in% c("blank", "not_stained")) %>%
  filter(! is.na(subclone)) %>%
  select(sample, day, `Record Time`, corrected_count) %>%
  group_by(sample, day) %>%
  mutate(replicate = row_number()) %>%
  mutate(corrected_count = log2(corrected_count)) %>%
  ungroup() %>%
  pivot_wider(
    names_from = day,
    values_from = c(`Record Time`, corrected_count),
    names_sep = "_day"
  ) %>%
  mutate(
    elapsed_hours = as.numeric(
      difftime(`Record Time_day2`, `Record Time_day1`, units = "hours")
    ),
    log2_growth = corrected_count_day2 - corrected_count_day1,
    doubling_time_hours = elapsed_hours / log2_growth,
    doubling_time_days = doubling_time_hours / 24
  ) %>%
  left_join(
    wells_metadata %>%
      select(sample, subclone, passage, aneu_group) %>%
      distinct(),
    by = "sample"
  ) 

doubling_time %>%
  filter(subclone != "mother clone" & passage == 2) %>%
  arrange(aneu_group, subclone) %>%
  mutate(subclone = factor(subclone, levels = unique(subclone))) %>%
  ggplot(aes(x = subclone, y = doubling_time_hours, group = sample, fill = aneu_group))+
  geom_boxplot()+
  geom_jitter()+
  scale_fill_manual(values = aneu_colors)



