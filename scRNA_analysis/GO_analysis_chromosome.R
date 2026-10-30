#Script wrote to do GO analysis of marker genes per chromosome
#load libraries
library(tidyverse)
source("git_modules/GO_tools/GO_tools.R")

#set paths
out_dir <- "../outputs/"

#get the marker genes
marker_genes <- readxl::read_xlsx("../outputs/scRNA/default_seurat_processing/cluster_markers.xlsx")
GO_terms <- read.delim("../inputs/reference_genome/LdCL_GO_terms_updated.tsv") %>%
  select(gene_id, simplified_term)

#get metadata for GO terms
GO_data <- get_GO_data(GO_terms = GO_terms$simplified_term)

#combine the data into a single data frame
GO_terms <- GO_terms %>%
  left_join(
    marker_genes %>%
      rename(gene_id = gene_id.x) %>%
      select(gene_id, cluster) %>%
      group_by(gene_id) %>%
      mutate(cluster = paste0(cluster, collapse = ",")) %>%
      distinct()
    ,
    by = "gene_id"
) %>%
  filter(! is.na(cluster)) %>%
  mutate(chromosome = as.numeric(gsub("^LdCL_([0-9]{2}).*$", "\\1", gene_id))) %>%
  left_join(
    GO_data %>%
      rename(simplified_term = input_term) %>%
      select(simplified_term, aspect, name)
    , by = "simplified_term"
  ) %>%
  left_join(
    marker_genes %>%
      rename(gene_id = gene_id.x, gene_name = cellranger_name) %>%
      select(gene_id, gene_name) %>%
      distinct()
    , by = "gene_id"
  )

#perform GO analysis per chromosome
GO_results <- list()
for(Chromosome in unique(GO_terms$chromosome)){
  GO_results[[Chromosome]] <- list()
    for(Aspect in unique(GO_terms$aspect)){
    set_terms <- GO_terms %>%
      filter(aspect == Aspect & chromosome == Chromosome) %>%
      pull("simplified_term")
    
    if(length(set_terms) == 0){
      next
    }else{
      GO_results[[Chromosome]][[Aspect]] <- GO_enrich(set_terms = set_terms, ref_terms = GO_terms$simplified_term, test_type = "enrichment") 
    }
    }
  GO_results[[Chromosome]] <- bind_rows(GO_results[[Chromosome]], .id = "aspect")
}
GO_results <- bind_rows(GO_results, .id = "chromosome")


#add metadata to GO terms
GO_data <- get_GO_data(GO_results$term)

GO_results <- GO_results %>%
  left_join(
    GO_data %>%
      rename(term = input_term, GO_name = name) %>%
      select(term, GO_name)
    , by = "term"
  )

#visualize it
plot_list <- list()
scale_labels <- format(10^-c(0:4), scientific = FALSE)
scale_labels[5] <- paste0("< ", scale_labels[5])
for(type in unique(GO_results$aspect)){
  plot_list[[type]] <- GO_results %>%
    filter(log2enrichment > 0) %>%
    filter(aspect == type)%>%
    mutate(aspect = gsub("_", " ", aspect)) %>%
    group_by(aspect) %>%
    mutate(GO_name = paste0(GO_name, ": (", term, ")")) %>%
    mutate(chromosome = factor(as.numeric(chromosome))) %>%
    mutate(GO_name = fct_reorder2(GO_name, dplyr::desc(pvalue), chromosome)) %>%
    mutate(pvalue = -log10(pvalue)) %>%
    mutate(pvalue = ifelse(pvalue > 4, 4, pvalue)) %>%
    ggplot(aes(x = chromosome, y = GO_name, fill = pvalue))+
    geom_tile()+
    labs(fill = "p-value", x = "cell cluster")+
    guides(fill = guide_colorbar(reverse = TRUE))+
    facet_grid(rows = vars(aspect), scales = "free", space = "free")+
    scale_fill_gradientn(colors = c("white", "purple"), limits = c(0, 4), breaks = c(0:4), labels = scale_labels)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    theme(panel.background = element_blank(), axis.title.y = element_blank())
}
panel_props <- unlist(lapply(plot_list, function(x)length(unique(x$data$term))))
panel_props <- sqrt(panel_props)
panel_props <- panel_props/sum(panel_props)

final_plot <- cowplot::plot_grid(plotlist = plot_list, ncol = 1, align = "vh", rel_heights = panel_props)


out_dir <- paste0(out_dir, "/Figures/")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

plot_scale <- 1.8
ggsave(paste0(out_dir, "/chromosomes_GO.pdf"), plot = final_plot, width = 6 * plot_scale, height = 10 * plot_scale)

#check which genes reger to GO tetrapyrrole byosynthetic process
GOI <- "GO:0033014"
GO_terms %>%
  filter(chromosome == 20 & simplified_term %in% GOI)