#libraries####
library(flowCore)
library(flowWorkspace)
library(openCyto)
library(tidyverse)
library(scales)
library(mclust)
library(ggridges)


#Analysis####
##load dataset####
fs <- read.flowSet(path = "inputs/growth_curves/FCS_files/BPK081_subclones_growth_curve_day14_200 ul/", pattern = "*.fcs", transformation = FALSE, truncate_max_range = FALSE)
##remove unecessary channels####
keep <- c("FSC-A","FSC-H","SSC-A","SSC-H","Time","FL2-A","FL10-A")
fs <- fs[,keep]

#transform fluorescent channels to log (keep FSC, SSC, and Time linear)####
channels <- colnames(fs)
fluo_channels <- channels[!grepl("^(FSC|SSC|Time)", channels, ignore.case = TRUE)] 
fluo_channels <- fluo_channels[grepl("-A", fluo_channels)] #keeps only Area parameter
#trans <- estimateLogicle(fs[[1]], channels = fluo_channels, m = 6)
trans <- transformList(fluo_channels, logTransform(transformation = "log10"))
fs <- transform(fs, trans)

#create a joint dataframe####
flow_df <- lapply(1:length(fs), function(i){
  df <- as.data.frame(flowCore::exprs(fs[[i]]))
  df$sample <- sampleNames(fs)[i]
  return(df)
})

flow_df <- bind_rows(flow_df)

#start gating####
##VybrantDyeCycle Vs SSC
#now visualize the Vybrant Dye Cycle signal

to_plot <- flow_df %>%
  filter(grepl("Well-A|Well-G", sample)) %>%
  mutate(`SSC-A` = ifelse(`SSC-A` > 2e6, 2e6, `SSC-A`)) %>% #cap the y-axis
  mutate(`FL2-A` = ifelse(`FL2-A` > 2e6, 2e6, `FL2-A`))

p <- to_plot %>%
  filter(grepl("Well-G3", sample)) %>%
  ggplot(aes(x = `FL2-A`, y = `SSC-A`))+
  geom_point(size = 0.05)+
  geom_hex(bins = 128)+
  scale_x_continuous(labels = log10_labels)+
  scale_y_continuous(labels = log10_labels)+
  scale_fill_viridis_c(option = "turbo")+
  guides(fill = "none")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  facet_wrap(vars(sample))

gate <- gate_from_ggplot(p, gate_name = "has_DNA")

p + geom_path(data = as.data.frame(gate@boundaries))

#visualize flow rate####
bins <- 500
flow_rate_plot <- flow_df %>%
  select(sample, Time) %>%
  #filter(grepl("G11", sample)) %>%
  group_by(sample) %>%
  mutate(bin = cut_interval(Time, n = bins, labels = FALSE)) %>%
  group_by(sample, bin) %>%
  mutate(bin_start = min(Time), bin_end = max(Time), bin_mid = median(c(bin_start, bin_end))) %>%
  summarise(bin_count = n(), bin_start = bin_start[1], bin_end = bin_end[1], bin_mid = bin_mid[1]) %>%
  group_by(sample) %>%
  mutate(thres_count = ifelse(bin_count <= 1, NA, bin_count)) %>%
  mutate(threshold = quantile(thres_count, 0.25, na.rm = TRUE)) %>%
  mutate(to_keep = bin_count > threshold) %>%
  ggplot(aes(x = bin_mid, y = bin_count, fill = to_keep))+
  geom_col()+
  geom_hline(aes(yintercept = threshold))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  facet_wrap(vars(sample), ncol = 1, scales = "free_y", strip.position = "left")+
  theme(
    panel.spacing.y = unit(-0.4, "lines"),  # overlap panels (ridge feel)
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0),
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    panel.border = element_blank()
  )

ggsave("outputs/growth_curves/QC/flow_rate.pdf", flow_rate_plot, width = 9, height = 30, limitsize = FALSE, create.dir = TRUE)


  
#now create a gate to separate events containing DNA (vybrant positive) from non-sepecific events
gs <- GatingSet(fs)