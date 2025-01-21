# usr/bin/Rscript
# %%
library(tidyverse)
library(dplyr)

known_loci_path <- c("/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_known_loci/datasrc/cteph_reported_loci.xlsx")
known_loci <- readxl::read_excel(known_loci_path)

graph_path <- c("/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_graph/graph/top10_gsea")
data_path <- c("/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_graph/data")

load("/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_graph/data/merged_pqtl_additive_apt.rda")

pqtl_list <- list()
for (snpid in ls()) {
  if (startsWith(snpid, "chr")) {
    pqtl_list[[snpid]] <- get(snpid)
  }
}

names(pqtl_list)

#%%
env_dfs <- new.env()

for (snpid in names(pqtl_list)) {
  ex_data <- pqtl_list[[snpid]]
  df <- as.data.frame(ex_data)
  env_dfs[[snpid]] <- df
}

# %%
snpids <- names(env_dfs)

# snpids <- c('chr1:169514323:T:C', 'chr9:133261703:A:G')
# snpids <- c('chr4:186286227:C:T', 'chr19:10631494:A:G', 'chr20:23048087:G:A')
# snpids <- c('chr4:154586438:T:C', 'chr4:154590745:T:C', 'chr20:34984375:G:A')
# snpids <- c('chr4:154599778:G:A')

get_protein_lists <- function(env_dfs, snpid) {
  pqtl_df <- env_dfs[[snpid]]
  
  case_df <- pqtl_df %>% dplyr::select(case_additive_beta, UniProt)
  control_df <- pqtl_df %>% dplyr::select(control_additive_beta, UniProt)
  
  proteinList_case <- case_df$case_additive_beta
  names(proteinList_case) <- case_df$UniProt
  proteinList_case <- sort(proteinList_case, decreasing = TRUE)
  
  proteinList_control <- control_df$control_additive_beta
  names(proteinList_control) <- control_df$UniProt
  proteinList_control <- sort(proteinList_control, decreasing = TRUE)
  
  return(list(case = proteinList_case, control = proteinList_control))
}

protein_lists <- list()
for (snpid in snpids) {
  protein_lists[[snpid]] <- get_protein_lists(env_dfs, snpid)
}

# %%
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(parallel)

run_gsea_kegg <- function(protein_lists, snpid, group) {
  set.seed(1234)
  gsea_result <- gseKEGG(
    protein_lists[[snpid]][[group]],
    organism = 'hsa',
    keyType = 'uniprot',
    exponent = 1,
    minGSSize = 10, 
    maxGSSize = 500,
    eps = 0,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    verbose = TRUE, 
    use_internal_data = FALSE, 
    seed = TRUE, 
    by = 'fgsea',
    nPermSimple = 100000  
  )
  return(gsea_result)
}

gsea_kegg_results <- mclapply(snpids, function(snpid) {
  list(
    case = run_gsea_kegg(protein_lists, snpid, "case"),
    control = run_gsea_kegg(protein_lists, snpid, "control")
  )
}, mc.cores = detectCores())

names(gsea_kegg_results) <- snpids

# %%
saveRDS(gsea_kegg_results, file = file.path(data_path, "gsea_kegg_results.rds"))

# %%
library(cowplot)
library(patchwork)
library(ggplot2)
library(stringr)
library(forcats)
library(DOSE)
library(grid)
library(gridExtra)


create_final_plot <- function(snpid_focus, known_loci, gsea_kegg_results) {

  library(extrafont)
  font_import()

  # Set global plot style
  theme_set(theme_gray(base_family = "DejaVu Sans", base_size = 12))

  gene_focus <- known_loci$GENE[known_loci$ID %in% snpid_focus]

  case_result <- gsea_kegg_results[[snpid_focus]]$case
  control_result <- gsea_kegg_results[[snpid_focus]]$control

  nes_range <- range(c(case_result$NES, control_result$NES))

  # create an empty plot
  empty_plot <- ggplot() + 
    geom_blank() + 
    theme_void() + 
    ggtitle("")

  # plot case result firstly
  if (nrow(case_result) > 0) {
    p1 <- ggplot(case_result, showCategory = 10, aes(NES, fct_reorder(Description, NES))) + 
      geom_segment(aes(xend=0, yend = Description)) +
      geom_point(aes(color=p.adjust, size = Count)) +
      scale_color_gradientn(colours=c("#f7ca64", "#46bac2", "#7e62a3"), trans = "log10", guide=guide_colorbar(reverse=TRUE, order=1)) +
      scale_size_continuous(range=c(3, 10)) +
      theme_dose(12) +
      xlab("NES") +
      ylab("") +
      labs(NULL) +
      ggtitle("CTEPH") + 
      theme(
        text = element_text(face = "bold", size = 18), 
        plot.margin = margin(1, 1, 1, 1, "cm"), 
        axis.text = element_text(size = 20), 
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20), 
        axis.text.y = element_text(size = 17),
        axis.text.x = element_text(size = 15)
      ) + 
      scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
      coord_cartesian(xlim = nes_range + c(-0.2, 0.2))
  } else {
    p1 <- empty_plot
  }

  # plot control result secondly
  if (nrow(control_result) > 0) {
    p2 <- ggplot(control_result, showCategory = 10, aes(NES, fct_reorder(Description, NES))) + 
      geom_segment(aes(xend=0, yend = Description)) +
      geom_point(aes(color=p.adjust, size = Count)) +
      scale_color_gradientn(colours=c("#f7ca64", "#46bac2", "#7e62a3"), trans = "log10", guide=guide_colorbar(reverse=TRUE, order=1)) +
      scale_size_continuous(range=c(3, 10)) +
      theme_dose(12) +
      xlab("NES") +
      ylab("") +
      labs(NULL) +
      ggtitle("Health") + 
      theme(
        text = element_text(face = "bold", size = 18), 
        plot.margin = margin(1, 1, 1, 1, "cm"), 
        axis.text = element_text(size = 20), 
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 17), 
        axis.text.x = element_text(size = 15)
      ) + 
      scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
      coord_cartesian(xlim = nes_range + c(-0.2, 0.2))
  } else {
    p2 <- empty_plot
  }

  # Create a title with background and border
  title <- ggdraw() + 
    draw_label(
      bquote(atop(italic(.(gene_focus)) ~ "--" ~ .(snpid_focus))),
      fontface = 'bold',
      fontfamily = "DejaVu Sans",
      size = 25, 
      vjust = 0.8
    ) +
    theme(
      plot.background = element_rect(fill = "lightblue", color = "lightblue", size = 1),
      plot.margin = margin(10, 10, 10, 10)
    )

  # arrange two plots with centered title
  final_plot <- plot_grid(p1, p2, ncol = 2, align = 'v', rel_widths = c(1, 1))
  final_plot <- plot_grid(title, final_plot, ncol = 1, rel_heights = c(0.1, 1))
  final_plot <- final_plot + theme(text = element_text(family = "DejaVu Sans"))
  
  return(final_plot)
}

# %%
library(extrafont)

options(ask = FALSE)

# Check if fonts are already imported
if (!"DejaVu Sans" %in% fonts()) {
  font_import(pattern = "DejaVu Sans", prompt = FALSE)
  loadfonts(device = "pdf")
}

snpids <- names(gsea_kegg_results)
for (snpid in snpids) {
  gene_name <- known_loci$GENE[known_loci$ID %in% snpid]
  final_plot <- create_final_plot(snpid, known_loci, gsea_kegg_results)
  ggsave(filename = file.path(graph_path, paste0(gene_name, "_", snpid, ".pdf")), plot = final_plot, device = "pdf", width = 20, height = 10)
}


# %%

get_common_gsea_results <- function(case_result, control_result) {
  # get common IDs
  case_ids <- case_result@result$ID
  control_ids <- control_result@result$ID
  common_ids <- intersect(case_ids, control_ids)
  # subset results
  case_result_sub <- case_result
  control_result_sub <- control_result
  
  case_result_sub@result <- case_result@result[case_result@result$ID %in% common_ids, ]
  control_result_sub@result <- control_result@result[control_result@result$ID %in% common_ids, ]
  
  return(list(case_result_sub = case_result_sub, control_result_sub = control_result_sub))
}


common_results <- get_common_gsea_results(case_result, control_result)
case_result_sub <- common_results$case_result_sub
control_result_sub <- common_results$control_result_sub

cat("Number of entries in case_result_sub:", nrow(case_result_sub@result), "\n")
cat("Number of entries in control_result_sub:", nrow(control_result_sub@result), "\n")


# %%
ids <- case_result_sub@result$ID

# Reorder control_result_sub to match the order of case_result_sub

ridgeplot(case_result_sub, core_enrichment = TRUE)
ridgeplot(control_result_sub, core_enrichment = TRUE)
