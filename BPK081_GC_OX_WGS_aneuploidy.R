#libraries
library(tidyverse)
library(readr)
library(mgcv)
library(ggalign)

#set output directory
out_dir <- "outputs/WGS/"

#inport bed files
file_names <- list.files("inputs/WGS/", pattern = "*counts.bed", recursive = TRUE)
files <- paste0("inputs/WGS/", file_names)
files <- lapply(files, read_delim)
sample_names <- gsub(".counts|.bed", "", file_names)
names(files) <- sample_names

#now we append the gc_content and compensated it
bin_meta <- read_delim("inputs/WGS/E4_G3_pTB007/count_files/bin_gc_and_mappability_1000.tsv")
files <- lapply(files, function(x){
  colnames(x) <- c("chrom", "start", "end", "read_count")
  x <- x %>%
    left_join(bin_meta, by = c("chrom", "start", "end")) %>%
    rename(chromosome = chrom,
           gc_content = `5_pct_gc`) %>%
    group_by(chromosome) %>%
    mutate(is_empty = read_count < 10) %>%
    mutate(is_outlier = read_count < quantile(read_count, 0.05) | read_count > quantile(read_count, 0.95)) %>%
    rowwise() %>%
    mutate(to_remove = any(is_empty, is_outlier, mappability < 1)) %>%
    ungroup()
  
  #to calculate gc bias we need to exclude aneuploid chromosomes. So we first calculate roughly which are the disomic chromosomes
  base_chromo <- x %>% 
    filter(!to_remove) %>%
    mutate(norm_read_count = read_count/median(read_count)) %>%
    group_by(chromosome) %>%
    summarise(mean_norm_read_count = median(norm_read_count)) %>%
    mutate(is_aneuploid = mean_norm_read_count < 0.9 | mean_norm_read_count > 1.1) %>%
    filter(!is_aneuploid) %>%
    pull(chromosome)
  
  #now we use gam modeling to determine the relationship between GC content and read count
  to_model <- x %>%
    filter(!to_remove) %>%
    filter(chromosome %in% base_chromo) %>%
    filter(gc_content < quantile(gc_content, 0.99) & gc_content > quantile(gc_content, 0.01)) %>%
    select(chromosome, read_count, gc_content)
    
  model <- gam(read_count ~ s(gc_content), data = to_model)
  
  #finally, we apply the model to correct the counts
  x <- x %>%
    mutate(gc_predicted = predict(model, as.data.frame(gc_content))) %>%
    mutate(gc_correction_factor = median(gc_predicted) / gc_predicted) %>%
    mutate(read_count_corrected = read_count * gc_correction_factor) %>%
    mutate(read_count_corrected = ifelse(read_count_corrected < 0, 0, read_count_corrected))
  
  return(x)
})
rm(bin_meta)

#bind the list as a single data frame
to_plot <- bind_rows(files, .id = "sample") %>%
  mutate(chromosome = as.integer(factor(chromosome))) %>%
  mutate(chromosome = ifelse(chromosome < 10, paste0("Ld0", chromosome), paste0("Ld", chromosome))) %>%
  group_by(sample) %>%
  arrange(chromosome, start, .by_group = TRUE) %>%
  mutate(bin = row_number()) %>%
  group_by(sample) %>%
  mutate(norm_read_count = ifelse(to_remove, NA, read_count_corrected)) %>%
  mutate(norm_read_count = norm_read_count/median(norm_read_count[which(chromosome == "Ld36")], na.rm = TRUE)) %>% #since chromosome 36 doesn't change in somy here, we used it for normalization. 
  ungroup()

#start by visualizing the effect of gc on coverage
plot_list <- list()
plot_list[["gc_bias_before"]] <- to_plot %>%
  group_by(sample) %>%
  mutate(norm_read_count = read_count/median(read_count)) %>% #to include noisy bins
  ggplot(aes(x = gc_content, y = norm_read_count, color = mappability == 1))+
  geom_point(size = 0.1)+
  geom_smooth(data = ~ filter(.x, mappability == 1), color = "blue")+
  scale_x_continuous(limits = c(0.4, 0.7))+
  scale_y_continuous(limits = c(0, 3))+
  scale_color_manual(values = c("darkred", "black"))+
  facet_wrap(vars(sample), nrow = 1)+
  guides(color = "none")+
  labs(title = "evaluation of GC bias: before correction (red = low-mappability bin)", x = "GC %", y = "Normalized read count")

plot_list[["gc_bias_after"]] <- to_plot %>%
  ggplot(aes(x = gc_content, y = norm_read_count))+
  geom_point(size = 0.1)+
  geom_smooth(color = "blue")+
  scale_x_continuous(limits = c(0.4, 0.7))+
  scale_y_continuous(limits = c(0, 3))+
  facet_wrap(vars(sample), nrow = 1)+
  labs(title = "evaluation of GC bias: after correction and removal of noisy bins", x = "GC %", y = "Normalized read count (corrected)")

#calculate somies
somies <- to_plot %>%
  filter(!to_remove) %>%
  group_by(sample, chromosome) %>%
  summarise(somy = median(norm_read_count) * 2)

#plot somies
heat_col <- c("#001221", "#002342", "#002342", "#014175", "#035ba3", "#00c3ff", "#00ffee", "#33ff00", "#ccff00", "#fffa00","#ffa600", "#D73027", "#A50026", "#541b1b", "#4d0600")

plot_list[["somies"]] <- somies %>%
  ggplot(aes(x = sample, y = fct_rev(chromosome), fill = somy, label = round(somy, 1)))+
  geom_tile()+
  geom_text(color = "black", size = 3.2)+
  geom_text(color = "white", size = 3)+
  scale_fill_gradientn(colors = heat_col, limits = c(0, 8), breaks = c(0:8))+
  labs(y = "Chromosome", title = "bulk somies")

#visualize coverage
plot_list[["coverage"]] <- to_plot %>%
  filter(!to_remove) %>%
  group_by(sample, chromosome) %>%
  mutate(norm_read_count = norm_read_count * 2) %>%
  mutate(median_count = median(norm_read_count)) %>%
  ggplot(aes(x = start, y = norm_read_count))+
  geom_point(size = 0.1)+
  geom_line(aes(y = median_count), color = "darkred")+
  facet_grid(rows = vars(sample), cols = vars(chromosome), scale = "free_x", switch = "x")+
  scale_x_continuous(breaks = NULL, labels = NULL, name = NULL)+
  scale_y_continuous(breaks = seq(0, 10, length = 6))+
  labs(y = "Normalized read count * 2", title = "coverage per 1kb bin")+
  theme_void()+
  theme(
    strip.clip = "off", 
    strip.text.x = element_text(angle = 45, hjust = 1),
    strip.text.y = element_text(angle = 270),
    axis.text.y = element_text(),
    axis.title.y = element_text(angle = 90),
  )

#build figure panel
left <-  ggalign::align_plots(plot_list[["gc_bias_before"]], plot_list[["gc_bias_after"]], ncol = 1) + layout_tags(tags = NULL)
up <- ggalign::align_plots(left, plot_list[["somies"]], nrow = 1, widths = c(0.7, 0.3))
final <- ggalign::align_plots(up, plot_list[["coverage"]], ncol = 1, heights = c(0.3, 0.7)) + layout_tags(tags = "A")

#scale text
plot_list <- lapply(plot_list, function(x)x+theme(text = element_text(size = 26)))

#save the final figure panel####
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
final_file <- paste0(out_dir, "/aneuploidy.pdf")
plot_scale <- 3
ggsave(final_file, plot = final, width = 8.27 * plot_scale, height = 7 * plot_scale)
ggsave(gsub(".pdf", ".png", final_file), plot = final, width = 8.27 * plot_scale, height = 7 * plot_scale, dpi = 600)

#open it
if (Sys.info()["sysname"] == "Darwin") {
  system2("open", args = final_file, wait = FALSE)
} else if (Sys.info()["sysname"] == "Linux") {
  system2("xdg-open", args = final_file, wait = FALSE)
}

