######################################################################
# Genome-Resolved Metagenomic Analysis (Country-Level)
# Alpha & Beta Diversity
######################################################################

set.seed(42)

# ---- Packages ----
suppressPackageStartupMessages({
  library(vegan)
  library(ggplot2)
  library(dplyr)
  library(FSA)         # Dunn test
  library(tidyr)
  library(tibble)
  library(parallel)
  library(RColorBrewer)
  library(ComplexHeatmap)
  library(pheatmap)
})

######################################################################
# 1. ALPHA DIVERSITY (raw TPM)
######################################################################

tpm <- read.csv("TPM_matrix.csv", row.names = 1, check.names = FALSE)
meta <- read.csv("Metadata.csv")

if (nrow(tpm) > ncol(tpm)) tpm <- t(tpm)
tpm[is.na(tpm)] <- 0

stopifnot("Sample_ID" %in% colnames(meta))
rownames(meta) <- meta$Sample_ID
meta <- meta[rownames(tpm), , drop = FALSE]

alpha_df <- data.frame(
  Sample_ID = rownames(tpm),
  Shannon   = diversity(tpm, "shannon"),
  Simpson   = diversity(tpm, "simpson"),
  Richness  = specnumber(tpm)
)

thresholds <- c(0.1, 1, 10)
for (thr in thresholds) {
  alpha_df[[paste0("Richness_TPM≥", thr)]] <- rowSums(tpm >= thr)
}
alpha_df <- left_join(alpha_df, meta, by = "Sample_ID")

# ---- Define color palette for countries ----
countries <- sort(unique(alpha_df$Country))
country_colors <- setNames(brewer.pal(min(length(countries), 8), "Set2"), countries)

# ---- Plotting Function ----
plot_alpha <- function(df, yvar, ylab, file) {
  p <- ggplot(df, aes(x = Country, y = .data[[yvar]], fill = Country)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
    scale_fill_manual(values = country_colors, drop = FALSE) +
    theme_minimal(base_size = 11) +
    labs(x = "", y = ylab) +
    theme(
      legend.position   = "none",
      axis.text.x       = element_text(angle = 45, hjust = 1, color = "black", size = 7),
      axis.text.y       = element_text(color = "black", size = 7),
      axis.title.y      = element_text(size = 9, face = "bold"),
      panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.ticks        = element_line(color = "black", linewidth = 0.3)
    )
  ggsave(file, p, width = 4, height = 3, dpi = 600)
}

# ---- Generate Plots ----
plot_alpha(alpha_df, "Shannon", "Shannon Index", "Alpha_Shannon_Country.tiff")
plot_alpha(alpha_df, "Simpson", "Simpson Index", "Alpha_Simpson_Country.tiff")
plot_alpha(alpha_df, "Richness", "Richness", "Alpha_Richness_Country.tiff")
for (thr in thresholds) {
  yvar <- paste0("Richness_TPM≥", thr)
  plot_alpha(alpha_df, yvar, paste0("Richness (TPM ≥ ", thr, ")"),
             paste0("Alpha_", yvar, "_Country.tiff"))
}

write.csv(alpha_df, "AlphaDiversity_Country.csv")

######################################################################
# 2. ALPHA DIVERSITY STATISTICS (Kruskal–Wallis + Dunn)
######################################################################

stat_list <- list()

kw_shannon  <- kruskal.test(Shannon ~ Country, data = alpha_df)
kw_simpson  <- kruskal.test(Simpson ~ Country, data = alpha_df)
kw_richness <- kruskal.test(Richness ~ Country, data = alpha_df)

kw_df <- tibble(
  Test        = "Kruskal-Wallis",
  Metric      = c("Shannon", "Simpson", "Richness"),
  p_value     = c(kw_shannon$p.value, kw_simpson$p.value, kw_richness$p.value)
)
stat_list[["KW"]] <- kw_df

# ---- Dunn’s Post Hoc ----
dunn_shannon  <- dunnTest(Shannon ~ Country, data = alpha_df, method = "bh")$res %>% mutate(Metric = "Shannon")
dunn_simpson  <- dunnTest(Simpson ~ Country, data = alpha_df, method = "bh")$res %>% mutate(Metric = "Simpson")
dunn_richness <- dunnTest(Richness ~ Country, data = alpha_df, method = "bh")$res %>% mutate(Metric = "Richness")

dunn_all <- bind_rows(dunn_shannon, dunn_simpson, dunn_richness) %>%
  select(Metric, Comparison, Z, P.adj)

alpha_stats <- bind_rows(stat_list[["KW"]], dunn_all)
write.csv(alpha_stats, "AlphaDiversity_Stats_Country.csv", row.names = FALSE)

######################################################################
# 3. BETA DIVERSITY (CLR-corrected TPM)
######################################################################

species_clr <- read.csv("TPM_clr_batch_corrected_v1.csv", row.names = 1, check.names = FALSE)
metadata <- read.csv("Metadata.csv")

if (nrow(species_clr) > ncol(species_clr)) species_clr <- t(species_clr)
rownames(metadata) <- metadata$Sample_ID
metadata <- metadata[rownames(species_clr), , drop = FALSE]

# ---- PERMANOVA ----
dist_matrix <- vegdist(species_clr, method = "euclidean")
anosim_res  <- anosim(dist_matrix, metadata$Country, permutations = 999)
print(anosim_res)

NMDS <- metaMDS(species_clr, k = 2, distance = "euclidean", trymax = 50)
nmds_df <- as.data.frame(NMDS$points) %>%
  rownames_to_column("Sample_ID") %>%
  left_join(metadata, by="Sample_ID")

# ---- NMDS Plot ----

library(ggplot2)
library(ggrepel)
library(dplyr)

# Load NMDS coordinates and metadata (as before)
# nmds_df must contain: MDS1, MDS2, Country, Continent (optional)

# Define color palette (consistent with previous figures)
country_colors <- c(
  "Benin" = "#1b9e77",
  "Burkina Faso" = "#d95f02",
  "China" = "#7570b3",
  "Finland" = "#66a61e",
  "Singapore" = "#e6ab02",
  "South Africa" = "#a6761d",
  "United Kingdom" = "#e7298a",
  "USA" = "#666666"
)

# Add NMDS centroids for clarity
centroids <- nmds_df %>%
  group_by(Country) %>%
  summarise(MDS1 = mean(MDS1), MDS2 = mean(MDS2))

# Plot
p <- ggplot(nmds_df, aes(x = MDS1, y = MDS2, color = Country)) +
  geom_point(size = 2.5, alpha = 0.75) +
  stat_ellipse(aes(fill = Country), geom = "polygon", alpha = 0.12, color = NA) +
  geom_text_repel(
    data = centroids, aes(label = Country),
    size = 3.3, fontface = "bold", color = "black", max.overlaps = 10
  ) +
  scale_color_manual(values = country_colors) +
  scale_fill_manual(values = country_colors) +
  theme_minimal(base_size = 13) +
  labs(
    x = "NMDS1",
    y = "NMDS2",
    title = "Country-Level Beta Diversity (CLR–Euclidean)",
    subtitle = "Species-level abundance profiles across built-environment microbiomes"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 11, color = "black"),
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray30"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2)
  )

ggsave("NMDS_Country_Level_withLegend.tiff", p, width = 6.5, height = 5.5, dpi = 600)
p

# ---- PERMANOVA (country-level) ----
perm <- adonis2(dist_matrix ~ Country, data = metadata, permutations = 999)
print(perm)
write.csv(as.data.frame(perm), "PERMANOVA_Country_CLR_Euclidean.csv")

######################################################################
# 4. Save Centroids (for Country Visualization)
######################################################################
dispersion <- betadisper(dist_matrix, metadata$Country)
centroids_df <- data.frame(Country = rownames(scores(dispersion, display = "centroids")),
                           scores(dispersion, display = "centroids"))
write.csv(centroids_df, "Country_Centroids.csv", row.names = FALSE)
#######################################







######################################################################
# Sample Distribution by Country and Environment
######################################################################
library(ggplot2)
library(dplyr)
library(RColorBrewer)

meta <- read.csv("Metadata.csv")

# Summarize counts
meta_summary <- meta %>%
  group_by(Group, Country) %>%
  summarise(Sample_Count = n(), .groups = "drop")

# Ensure consistent ordering
env_order <- c("Hosp_env", "Hosp_sewage", "Ambulance", "Public_transp")
meta_summary$Group <- factor(meta_summary$Group, levels = env_order)

# Define palette for countries
countries <- sort(unique(meta_summary$Country))
country_colors <- setNames(brewer.pal(min(length(countries), 8), "Set2"), countries)

# Create horizontal barplot
p <- ggplot(meta_summary, aes(x = Sample_Count, y = Group, fill = Country)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = country_colors, name = "Country") +
  theme_minimal(base_size = 12) +
  labs(
    x = "Number of Samples",
    y = "Environment Type",
    title = "Sample Composition by Environment and Country"
  ) +
  theme(
    axis.text.y = element_text(size = 10, face = "bold", color = "black"),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.title.x = element_text(size = 11, face = "bold"),
    axis.title.y = element_text(size = 11, face = "bold"),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid.major.y = element_blank()
  )

ggsave("Environment_vs_Country_Sample_Distribution.tiff", p, width = 7, height = 4.5, dpi = 600)

p

#######################################################











