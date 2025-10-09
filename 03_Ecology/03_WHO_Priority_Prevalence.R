# ======================= WHO Priority Heatmaps (Publication-Ready) =======================
set.seed(42)

# ---- Packages ----
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr)
  library(ComplexHeatmap); library(circlize); library(grid)
})

# ---- Load data ----
tpm  <- read.csv("TPM_matrix.csv", row.names = 1, check.names = FALSE)
meta <- read.csv("Merged_metadata_all.csv", stringsAsFactors = FALSE)
stopifnot(all(c("Sample_ID","Group") %in% colnames(meta)))
rownames(meta) <- meta$Sample_ID

# Clean species names
rownames(tpm) <- gsub("_", " ", rownames(tpm))

# Align samples in TPM with metadata
common_samples <- intersect(colnames(tpm), meta$Sample_ID)
tpm  <- tpm[, common_samples, drop = FALSE]
meta <- meta[common_samples, , drop = FALSE]

# Replace NAs with 0
tpm[is.na(tpm)] <- 0

# ---- WHO lists ----
who_critical <- c("Acinetobacter baumannii","Pseudomonas aeruginosa",
                  "Klebsiella pneumoniae","Escherichia coli","Staphylococcus aureus")
who_high <- c("Enterococcus faecalis","Klebsiella oxytoca",
              "Stenotrophomonas maltophilia","Streptococcus pneumoniae")
who_other <- c("Proteus mirabilis","Enterococcus faecium","Neisseria gonorrhoeae",
               "Haemophilus influenzae","Helicobacter pylori","Clostridioides difficile",
               "Salmonella Typhi","Candida albicans","Candida parapsilosis",
               "Candida glabrata","Candida dubliniensis")
who_genus <- c("Enterobacter","Citrobacter","Serratia","Providencia","Morganella","Shigella","Campylobacter")

species_all <- rownames(tpm)
who_exact <- species_all[species_all %in% c(who_critical, who_high, who_other)]
who_from_genus <- unlist(lapply(who_genus, function(g)
  grep(paste0("^", g, "\\b"), species_all, value = TRUE, ignore.case = TRUE)))
who_list <- unique(c(who_exact, who_from_genus))
if (length(who_list) == 0) stop("No WHO species matched TPM row names.")

# ---- Map species to WHO categories ----
category_map <- setNames(rep("Other", length(who_list)), who_list)
category_map[names(category_map) %in% who_critical] <- "Critical"
category_map[names(category_map) %in% who_high]     <- "High"
category_map[names(category_map) %in% who_other]    <- "Other"

# Keep only WHO species
tpm_who <- tpm[rownames(tpm) %in% who_list, , drop = FALSE]
cat_order <- c("Critical","High","Other")
row_order <- order(factor(category_map[rownames(tpm_who)], levels = cat_order), rownames(tpm_who))
tpm_who   <- tpm_who[row_order, , drop = FALSE]
category_map <- category_map[rownames(tpm_who)]

# ---- Collapse to environments ----
env_order <- c("Ambulance","Hosp_env","Hosp_sewage","Public_transp")
env_order <- env_order[env_order %in% meta$Group]

log_tpm <- log10(tpm_who + 1)
mean_by_env <- do.call(cbind, lapply(env_order, function(grp) {
  idx <- which(meta$Group == grp)
  if (length(idx) == 0) return(rep(NA_real_, nrow(log_tpm)))
  rowMeans(log_tpm[, idx, drop = FALSE], na.rm = TRUE)
}))
colnames(mean_by_env) <- env_order
rownames(mean_by_env) <- rownames(log_tpm)

# ---- Save tables ----
write.csv(mean_by_env, "WHO_mean_logTPM_by_environment.csv")

# ---- Colors ----
env_cols <- c(Ambulance="#1f77b4", Hosp_env="#aec7e8",
              Hosp_sewage="#ff7f0e", Public_transp="#fdbf6f")
env_cols <- env_cols[colnames(mean_by_env)]
who_cols <- c(Critical="#d62728", High="#ff7f0e", Other="#7f7f7f")

# Define row splits and titles
row_split_vec <- factor(category_map, levels = cat_order)
split_titles  <- setNames(c("CRITICAL", "HIGH", "OTHER"), cat_order)

# Continuous color function
col_fun <- circlize::colorRamp2(c(0,1,2,3,4),
                                c("#f7fbff","#c6dbef","#6baed6","#2171b5","#08306b"))

# ---- Annotations ----
col_anno <- HeatmapAnnotation(
  Environment = colnames(mean_by_env),
  col = list(Environment = env_cols),
  annotation_name_gp = gpar(fontsize = 0),
  border = TRUE
)
row_anno <- rowAnnotation(
  WHO_Category = row_split_vec,
  col = list(WHO_Category = who_cols),
  annotation_name_gp = gpar(fontsize = 0),
  width = unit(4, "mm"), border = TRUE
)

# ---- Final Figure ----
tiff("WHO_heatmap_meanTPM_clean.tiff", width = 7, height = 5,
     units = "in", res = 600, compression = "lzw")
ht <- Heatmap(
  mean_by_env,
  name = "log10(TPM+1)",
  col = col_fun,
  top_annotation = col_anno,
  left_annotation = row_anno,
  row_split = row_split_vec,
  row_title = split_titles[levels(row_split_vec)],
  row_title_gp = gpar(fontsize = 8, fontface = "bold"),
  cluster_rows = FALSE, cluster_columns = FALSE,
  show_row_names = TRUE, row_names_gp = gpar(fontsize = 8),
  row_names_max_width = unit(20, "cm"),
  show_column_names = TRUE, column_names_gp = gpar(fontsize = 8, fontface = "bold"),
  column_names_rot = 90,
  border = TRUE,
  heatmap_legend_param = list(
    title_gp = gpar(fontface="bold"),
    at = 0:4, labels = 0:4,
    legend_direction = "horizontal",
    legend_width = unit(4, "cm")
  )
)
draw(ht, padding = unit(c(8,8,8,8), "mm"),
     heatmap_legend_side = "bottom", annotation_legend_side = "right")
dev.off()
####################################################
#Calculating the Prevalence

# ---- 2) Convert to presence/absence ----
pa_matrix <- tpm
pa_matrix[pa_matrix > 0] <- 1

pa_who <- pa_matrix[rownames(pa_matrix) %in% who_list, ]

# ---- Align metadata ----
meta <- meta[match(colnames(pa_who), meta$Sample_ID), ]


# ---- Calculate prevalence (%) per environment ----
prevalence_who <- pa_who %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Species") %>%
  pivot_longer(-Species, names_to = "Sample_ID", values_to = "Present") %>%
  left_join(meta[, c("Sample_ID", "Group")], by = "Sample_ID") %>%
  group_by(Group, Species) %>%
  summarise(Prevalence = mean(Present) * 100, .groups = "drop")

# ---- Classify species ----
prevalence_who <- prevalence_who %>%
  mutate(Category = case_when(
    Prevalence >= 80 ~ "Core",
    Prevalence >= 50 ~ "Secondary",
    TRUE ~ "Peripheral"
  ))

# ---- Count per category/environment ----
summary_counts <- prevalence_who %>%
  group_by(Group, Category) %>%
  summarise(n_species = n(), .groups = "drop")

# ---- Save outputs ----
write.csv(prevalence_who, "WHO_species_prevalence_classification.csv", row.names = FALSE)
write.csv(summary_counts, "WHO_species_category_counts.csv", row.names = FALSE)

#####################################################################################
# Load packages
library(igraph)
library(ggraph)
library(ggplot2)

# ---- 1. Load TPM data ----
tpm <- read.csv("TPM_matrix.csv", row.names = 1, check.names = FALSE)
tpm[is.na(tpm)] <- 0   # Replace NAs with zero

# Transpose if rows = species instead of samples
if(nrow(tpm) > ncol(tpm)) {
  tpm <- t(tpm)
}

# ---- 2. Define WHO pathogens ----
exact_matches <- c(
  "Escherichia coli", "Acinetobacter baumannii", "Mycobacterium tuberculosis", "Salmonella Typhii",
  "Enterococcus faecium", "Pseudomonas aeruginosa", "Neisseria gonorrhoeae",
  "Staphylococcus aureus", "Streptococcus pneumoniae",
  "Haemophilus influenza", "Helicobacter pylori", "Clostridioides difficile",
  "Klebsiella oxytoca", "Enterococcus faecalis", "Enterococcus avium",
  "Candida albicans", "Candida parapsilosis", "Candida glabrata", "Candida dubliniensis",
  "Proteus mirabilis", "Klebsiella pneumoniae", "Stenotrophomonas maltophilia"
)

group_matches <- list(
  "Shigella" = "Shigella",
  "Enterobacter" = "Enterobacter",
  "Citrobacter" = "Citrobacter",
  "Proteus" = "Proteus",
  "Serratia" = "Serratia",
  "Streptococcus" = "Streptococcus",
  "Morganella" = "Morganella",
  "Providencia" = "Providencia",
  "Campylobacter" = "Campylobacter"
)

species_names <- colnames(tpm)
all_priority_species <- species_names[
  species_names %in% exact_matches |
    sapply(species_names, function(x) {
      any(sapply(group_matches, function(pattern) grepl(pattern, x, ignore.case = TRUE)))
    })
]

# ---- 3. Filter species (≥5% prevalence) ----
filtered_matrix <- tpm[, colSums(tpm > 0) > (0.05 * nrow(tpm))]

# ---- 4. Correlation matrix ----
cor_matrix <- cor(filtered_matrix, method = "spearman")
diag(cor_matrix) <- 0


# ---- Build igraph as before ----
threshold <- 0.5
adj_matrix <- ifelse(abs(cor_matrix) > threshold, cor_matrix, 0)  # keep sign
diag(adj_matrix) <- 0
network <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE)

# ---- Node attributes ----
V(network)$type <- ifelse(V(network)$name %in% all_priority_species, "WHO Priority", "Other")
V(network)$degree <- degree(network)
V(network)$size <- ifelse(V(network)$type == "WHO Priority", 8, 3) + sqrt(V(network)$degree)
V(network)$label <- ifelse(V(network)$type == "WHO Priority", V(network)$name, NA)

# ---- Edge attributes ----
E(network)$sign <- ifelse(E(network)$weight > 0, "Positive", "Negative")


############################################################
# FIGURE 1: Full Co-occurrence Network
############################################################
# ---- Plot full network ----
tiff("WHO_Priority_Cooccurrence_Network_upgraded.tiff", width = 6, height = 6, units = "in", res = 600, compression = "lzw")
ggraph(network, layout = "fr") +
  geom_edge_link(aes(color = sign, alpha = abs(weight)), show.legend = TRUE) +
  scale_edge_color_manual(values = c("Positive" = "steelblue", "Negative" = "tomato")) +
  scale_edge_alpha(range = c(0.2, 0.8)) +
  geom_node_point(aes(size = size, color = type)) +
  scale_color_manual(values = c("WHO Priority" = "red", "Other" = "grey70")) +
  geom_node_text(aes(label = label), repel = TRUE, size = 3, fontface = "bold") +
  theme_void() +
  theme(legend.position = "right") +
  ggtitle("")
dev.off()

############################################################
# FIGURE 2: WHO-only Subnetwork
############################################################
# ---- Plot sub network ----
who_nodes <- V(network)$name[V(network)$name %in% all_priority_species]
who_network <- induced_subgraph(network, vids = who_nodes)

p2 <- ggraph(who_network, layout = "fr") +
  geom_edge_link(aes(width = weight), color = "gray50", alpha = 0.6) +
  geom_node_point(aes(size = degree), color = "#d62728", alpha = 0.9) +
  geom_node_text(aes(label = name), repel = TRUE, size = 4, color = "black") +
  scale_size(range = c(4, 14)) +
  theme_void() +
  labs(title = "")

# Save
tiff("WHO_Priority_Subnetwork.tiff", width = 6, height = 6, units = "in", res = 600, compression = "lzw")
print(p2)
dev.off()

##########################################################################################333

############################################################
# WHO Priority Pathogen Co-occurrence Network Statistics
############################################################



species_names <- colnames(tpm)
all_priority_species <- species_names[
  species_names %in% exact_matches |
    sapply(species_names, function(x) {
      any(sapply(group_matches, function(pattern) grepl(pattern, x, ignore.case = TRUE)))
    })
]

# ---- 3. Build network ----
filtered_matrix <- tpm[, colSums(tpm > 0) > (0.05 * nrow(tpm))]
cor_matrix <- cor(filtered_matrix, method = "spearman")
diag(cor_matrix) <- 0

threshold <- 0.5   # adjustable
adj_matrix <- ifelse(abs(cor_matrix) > threshold, abs(cor_matrix), 0)

network <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE)

# ---- 4. Compute centrality measures ----
stats <- data.frame(
  Species = V(network)$name,
  Degree = degree(network),
  Betweenness = betweenness(network, normalized = TRUE),
  Closeness = closeness(network, normalized = TRUE),
  Eigenvector = evcent(network)$vector,
  ClusteringCoeff = transitivity(network, type = "local", isolates = "zero")
)

# ---- 5. Subset WHO pathogens only ----
who_stats <- stats[stats$Species %in% all_priority_species, ]

# ---- 6. Save ----
write.csv(who_stats, "WHO_Priority_Network_Statistics.csv", row.names = FALSE)

# ---- 7. Optional: Rank by importance ----
top_who <- who_stats[order(-who_stats$Eigenvector), ]
print(top_who)


# ---- 3. Figure A: Barplot of top hubs (Degree + Eigenvector) ----
top_hubs <- who_stats %>% arrange(desc(Degree)) %>% slice(1:10)

ggplot(top_hubs, aes(x = reorder(Species, Degree), y = Degree, fill = Eigenvector)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient(low = "lightblue", high = "darkred") +
  theme_minimal(base_size = 12) +
  labs(title = "WHO Priority Pathogen Hubs",
       x = "Species", y = "Degree (Connections)",
       fill = "Eigenvector\nCentrality")
ggsave("WHO_Hubs_Barplot.tiff", width = 7, height = 5, dpi = 600)

# ---- 4. Figure B: Radar Plot (Risk Profiles) ----
# Normalize metrics 0–1
scaled <- as.data.frame(scale(who_stats[, -1]))
scaled$Species <- who_stats$Species

# Select a few key pathogens for visualization
radar_subset <- scaled %>% filter(Species %in% c("Escherichia coli",
                                                 "Klebsiella pneumoniae",
                                                 "Staphylococcus aureus",
                                                 "Acinetobacter baumannii",
                                                 "Pseudomonas aeruginosa"))

# Radar format: fmsb needs max/min rows
radar_data <- rbind(apply(radar_subset[, -ncol(radar_subset)], 2, max),
                    apply(radar_subset[, -ncol(radar_subset)], 2, min),
                    radar_subset[, -ncol(radar_subset)])
rownames(radar_data) <- c("Max", "Min", radar_subset$Species)

# Plot radar
tiff("WHO_Risk_Radar.tiff", width = 7, height = 7, units = "in", res = 600)
radarchart(radar_data,
           axistype = 1,
           pcol = rainbow(nrow(radar_data)-2),
           plwd = 3,
           cglcol = "grey", cglty = 1,
           axislabcol = "black", caxislabels = seq(0, 1, 0.2),
           vlcex = 0.8,
           title = "WHO Pathogen Risk Profiles")
legend("topright", legend = radar_subset$Species,
       col = rainbow(nrow(radar_data)-2), lty = 1, lwd = 3, bty = "n")
dev.off()

# ---- 5. Figure C: Community Detection ----
comm <- cluster_louvain(as.undirected(network))
V(network)$community <- membership(comm)

# Highlight WHO pathogens
V(network)$type <- ifelse(V(network)$name %in% all_priority_species, "WHO", "Other")

tiff("WHO_Modularity_Communities.tiff", width = 8, height = 7, units = "in", res = 600)
ggraph(network, layout = "fr") +
  geom_edge_link(aes(alpha = ..index..), color = "grey70") +
  geom_node_point(aes(color = as.factor(community), size = ifelse(type == "WHO", 6, 2))) +
  geom_node_text(aes(label = ifelse(type == "WHO", name, "")),
                 repel = TRUE, size = 3) +
  theme_void() +
  labs(title = "WHO Priority Pathogens in Community Modules",
       color = "Community")
dev.off()
################################################################################
