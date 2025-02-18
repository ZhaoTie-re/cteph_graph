# usr/bin/Rscript
# %%
library(tidyverse)
library(dplyr)
library(DOSE)

enrich_results_path <- c("/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_graph/data/gsea_kegg_results.rds")
enrich_results <- readRDS(enrich_results_path)

#  %%
snpid <- 'chr20:34984375:G:A'
case_kegg <- enrich_results[[snpid]]$case
control_kegg <- enrich_results[[snpid]]$control

case_kegg_df <- as.data.frame(case_kegg)
control_kegg_df <- as.data.frame(control_kegg)

case_kegg_df <- case_kegg_df %>%
    dplyr::select(ID, Description, setSize, NES, pvalue) %>%
    rename('case_Description' = Description, 'case_setSize' = setSize, 'case_NES' = NES, 'case_pvalue' = pvalue)

control_kegg_df <- control_kegg_df %>% 
    dplyr::select(ID, Description, setSize, NES, pvalue) %>%
    rename('control_Description' = Description, 'control_setSize' = setSize, 'control_NES' = NES, 'control_pvalue' = pvalue)

merged_df <- merge(case_kegg_df, control_kegg_df, by = "ID", all = TRUE)

merged_df <- merged_df %>%
    mutate(Description = ifelse(is.na(case_Description), control_Description, case_Description)) %>%
    select(-case_Description, -control_Description) %>%
    mutate(setSize = ifelse(is.na(case_setSize), control_setSize, case_setSize)) %>%
    select(-case_setSize, -control_setSize) %>%
    mutate(group = case_when(
        !is.na(case_NES) & !is.na(control_NES) ~ "group1",
        !is.na(case_NES) & is.na(control_NES) ~ "group2",
        is.na(case_NES) & !is.na(control_NES) ~ "group3"
    )) %>%
    arrange(group)

# write.csv(merged_df, file = "merged_df.csv", row.names = FALSE)

# %%
known_loci_path = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_known_loci/datasrc/cteph_reported_loci.xlsx" 
known_loci_df <- readxl::read_excel(known_loci_path)

gene <- known_loci_df %>% filter(ID == snpid) %>% pull(GENE)

# %%
output_path <- "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_graph/graph/gsea"

# %%
library(ggplot2)
library(ggh4x)
library(dplyr)
library(tidyr)

long_df <- merged_df %>%
    gather(key = "type", value = "NES", case_NES, control_NES)

# add order column to long_df for ggplot to use for ordering
long_df <- long_df %>%
    mutate(group = case_when(
        group == "group1" ~ "Common in CTEPH and Healthy",
        group == "group2" ~ "Unique in CTEPH",
        group == "group3" ~ "Unique in Healthy",
        TRUE ~ group
    )) %>%
    group_by(group) %>%
    mutate(order = case_when(
        group == "Common in CTEPH and Healthy" & type == "case_NES" ~ dense_rank(desc(abs(NES))),
        group == "Unique in CTEPH" & type == "case_NES" ~ dense_rank(desc(abs(NES))),
        group == "Unique in Healthy" & type == "control_NES" ~ dense_rank(desc(abs(NES))),
        group == "Unique in Healthy" & type == "case_NES" ~ dense_rank(desc(abs(NES)))
    )) %>%
    ungroup() %>%
    group_by(group, ID) %>%
    mutate(order = case_when(
        group == "Common in CTEPH and Healthy" & type == "control_NES" ~ order[type == "case_NES"],
        group == "Unique in CTEPH" & type == "control_NES" ~ order[type == "case_NES"],
        group == "Unique in Healthy" & type == "case_NES" ~ order[type == "control_NES"],
        TRUE ~ order
    )) %>%
    ungroup()



# plot lollipop chart with fixed facet size
ggplot(long_df, aes(x = reorder(Description, order), y = NES, group = type)) +
    geom_point(aes(size = setSize, color = type), alpha = 0.8) +
    geom_point(aes(size = setSize), shape = 21, color = "black", fill = NA) +
    facet_wrap(~ group, scales = "free_x", ncol = 1) + 
    coord_cartesian(ylim = c(min(long_df$NES, na.rm = TRUE), max(long_df$NES, na.rm = TRUE))) +  
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 25, face = "bold"),
          axis.text.y = element_text(size = 30, face = "bold"),
          axis.title.y = element_text(size = 45, face = "bold", margin = margin(r = 200)), 
          panel.grid.major = element_line(color = "grey80"),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "lightblue", color = "black"),
          strip.text = element_text(face = "bold", size = 22),
          axis.title.x = element_blank(),
          plot.title = element_text(size = 40, face = "bold"), 
          legend.title = element_text(size = 18, face = "bold"),
          legend.text = element_text(size = 16, face = "bold"),
          plot.margin = unit(c(1, 5, 1, 15), "cm"),  
          panel.spacing = unit(2, "lines"),  
          strip.text.x = element_text(size = 25, face = "bold"),
          legend.justification = "center") + 
    labs(title = bquote(bolditalic(.(gene)) ~ bold("--") ~ bold(.(snpid))),
         y = "NES",
         size = "Set Size",
         color = "Group") +
    scale_size_continuous(range = c(5, 15)) +  # increase the range of point sizes
    guides(color = guide_legend(override.aes = list(size = 10))) +  # increase legend point size
    scale_color_manual(labels = c("CTEPH", "Healthy"), values = c("case_NES" = "blue", "control_NES" = "green"), name = "Group")


# %%
fold <- length(unique(long_df$group))
ggsave(paste0(output_path, "/gsea_", gene, "_", snpid, ".pdf"), width = 70, height = 10*fold, limitsize = FALSE)
