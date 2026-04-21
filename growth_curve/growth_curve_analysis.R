#libraries
library(tidyverse)
library(readODS)
library(ggalign)
#supporting functions
#parameters
inputs_dir <- "inputs/growth_curves/statistics/"
wells_metadata <- "inputs/growth_curves/wells_metadata.ods"
cols_to_plot <- c("Has DNA % Parent", 
                  "Is alive % Parent", 
                  "Singlets % Parent", 
                  "Singlets Median FSC-A", 
                  "Singlets Median B610-A", 
                  "Singlets GeoMean FSC-A",
                  "Singlets GeoMean B610-A", 
                  "Is alive SD Time",
                  "original_concentration(par/ml)")

#load wells metadata
wells_metadata <- read_ods(wells_metadata, sheet = 1)

#load files
file_names <- list.files(inputs_dir, pattern = "*.csv")
files <- paste0(inputs_dir, file_names)
files <- lapply(files, read.csv, skip = 2, header = FALSE) %>%
  lapply(
    function(x){
      colnames(x) <- x[1,] 
      x <- x[-1,]
      return(x)}
    )
names(files) <- file_names

#format dataframes
files <- lapply(files, function(x){
  #fix colnames
  colnames(x) <- gsub(":", "", colnames(x))
  
  #convert `####` to `NA`
  x[x == "####"] <- NA
  
  #convert characters to their corresponding datatypes
  x <- type.convert(x, as.is = TRUE, dec = ".")
  
  #convert percentages to numeric
  x_percent <- x[grep("%", x)]
  x_percent <- apply(x_percent, 2, function(y){
    y <- as.numeric(gsub("%", "", y))
    y <- y/100
  })
  x[grep("%", x)] <- x_percent
  rm(x_percent)
  

  
  #convert date to proper posix
  x[["Record Time"]] <- as.POSIXct(x[["Record Time"]],
                                   format = "%Y-%m-%d %H:%M:%S",
                                   tz = "UTC")
  
  
  #get row and column position of the wells
  x <- x %>%
    mutate(
      well = gsub("01-Well-","",`Tube Name`),
      row = str_extract(well, "[A-Z]"),
      column = as.numeric(str_extract(well, "\\d+$")))
  
  #convert rows and columns to factor
  x$row <- factor(x$row, LETTERS[1:8])
  x$column <- factor(x$column, 1:12)
  return(x)
})

#now append the list of objects into a single data frame
to_plot <- bind_rows(files, .id = "file")

#append the wells metadata
to_plot <- to_plot %>%
  left_join(wells_metadata, by = "well")

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
ggsave("outputs/growth_curves/plate_design.pdf", final_plot, width = 11.69*scale_factor, height = 8.27*scale_factor, create.dir = TRUE)

#plot the columns specified above
plot_list <- list()
for(col in cols_to_plot){
  plot <- to_plot %>%
    select(row, column, file, all_of(col)) %>%
    mutate(label = .data[[col]])
  
  #format label column
  if(grepl("%", col)){
    plot$label <- paste0(round(plot$label * 100, 1), "%")
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
  ggsave(paste0("outputs/growth_curves/", col, ".pdf"), plot + facet_wrap(vars(file)), width = 4+(1*ncols*scale_factor), height = 1*ncols*scale_factor, create.dir = TRUE)
}

final_plot <- ggalign::align_plots(!!!plot_list, ncol = 1)
scale_factor <- 2
ggsave("outputs/growth_curves/numeric_variables.pdf", final_plot, width = 11.69*scale_factor, height = 1*length(cols_to_plot)*scale_factor, create.dir = TRUE)

#now plot the counts per day
to_plot <- to_plot %>%
  mutate(day = as.Date(to_plot$`Record Time`, tz = "UTC")) %>%
  mutate(day = as.integer(difftime(day, min(day, na.rm = TRUE), units = "days")))

to_plot %>%  
  filter(!sample %in% c("blank", "not_stained")) %>%
  filter(! file %in% names(files)[c(4,7,10)]) %>%
  filter(! is.na(subclone)) %>%
  mutate(label = `original_concentration(par/ml)`) %>%
  mutate(label = format(label, scientific = TRUE, digits = 2)) %>%
  group_by(subclone, passage, day) %>%
  mutate(count = `original_concentration(par/ml)`) %>%
  mutate(
    med = median(count, na.rm = TRUE),
    mad0 = mad(count, constant = 1, na.rm = TRUE),
    outlier = if_else(mad0 > 0, abs(count - med) > 3 * mad0, FALSE)
  ) %>%
  #mutate(count = ifelse(outlier, NA, count)) %>%
  mutate(count = mean(count, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(!day %in% c(5)) %>%
  ggplot(aes(x = day, y = count, group = well, color = aneu_group))+
  geom_line()+
  facet_grid(rows = vars(subclone))
  
#quick check if using 200 ul instead of 100 ul affect counting
dil_diff <- to_plot %>%
  filter(day %in% c(4, 6)) %>%
  select(well, day, file, `original_concentration(par/ml)`) %>%
  group_by(day) %>%
  mutate(dilution = ifelse(grepl("200", file), "200ul", "100ul")) %>%
  select(well, day, dilution, `original_concentration(par/ml)`) %>%
  pivot_wider(names_from = dilution, values_from = `original_concentration(par/ml)`) %>%
  mutate(diff = abs(`200ul` - `100ul`)) %>%
  arrange(day, desc(diff))

dil_diff %>%
  ggplot(aes(x = `100ul`, y = `200ul`, label = well))+
  geom_text()+
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
  facet_grid(cols = vars(day))
