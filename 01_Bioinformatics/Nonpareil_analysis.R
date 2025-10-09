set.seed(45)
suppressPackageStartupMessages({
  stopifnot(requireNamespace("Nonpareil", quietly = TRUE))
  library(Nonpareil)
  library(dplyr)
  library(readr)
  library(tibble)
})

## folders on disk, also the labels you want to see on plots
projects <- c("Ambulance", "Hosp_env", "Hosp_sewage", "Public_transport")

## colors keyed by those exact labels
proj_cols <- c(
  "Ambulance"       = "#1f77b4",
  "Hosp_env"        = "#aec7e8",
  "Hosp_sewage"     = "#ff7f0e",
  "Public_transport"= "#fdbf6f"
)

## build samples from folders, attach label and color
samples <- bind_rows(lapply(projects, function(p) {
  files <- if (dir.exists(p)) list.files(p, pattern = "\\.npo$", full.names = TRUE) else character()
  if (!length(files)) return(NULL)
  tibble(
    File    = files,
    Name    = basename(files),
    Col     = unname(proj_cols[[p]]),
    Project = p     # legend label equals folder name
  )
}))

if (is.null(samples) || nrow(samples) == 0) {
  stop("No .npo files found in any of: ", paste(projects, collapse = ", "))
}

write.table(samples, "samples.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

## all-projects plot
png("nonpareil_all_projects.png", width = 6, height = 5, units = "in", res = 600)
nps <- Nonpareil.set(
  files     = samples$File,
  col       = samples$Col,
  labels    = samples$Name,
  plot      = TRUE,
  plot.opts = list(
    plot.observed   = TRUE,
    plot.model      = TRUE,
    plot.dispersion = FALSE,
    legend          = FALSE,
    lwd             = 1,
    main            = ""
  )
)
axis(1, lwd = 2.5, lwd.ticks = 1.2)
axis(2, lwd = 2.5, lwd.ticks = 1.2)
box(lwd = 2.5)

## legend, forced order matching your request, only those present are shown
leg_order <- projects
leg_labs  <- leg_order[leg_order %in% unique(samples$Project)]
legend("bottomright", inset = c(0.0045, 0.1),
       legend = leg_labs,
       col    = unname(proj_cols[leg_labs]),
       lwd    = 1,
       cex    = 0.7,
       bg     = "white",
       bty    = "n",
       title  = "")
dev.off()

## combined summary
sum_all <- as.data.frame(summary(nps)) %>%
  mutate(Name = samples$Name,
         File = samples$File,
         Project = samples$Project) %>%
  select(Project, Name, File, kappa, C, LR, LRstar, modelR, diversity)
write_csv(sum_all, "nonpareil_summary_all_projects.csv")

## per-project plots and summaries
dir.create("nonpareil_per_project", showWarnings = FALSE)
for (p in leg_labs) {  # iterate in the same order
  sub <- samples %>% filter(Project == p)
  out_png <- file.path("nonpareil_per_project", paste0(p, ".png"))
  png(out_png, width = 4, height = 4, units = "in", res = 600)
  Nonpareil.set(
    files     = sub$File,
    col       = sub$Col,
    labels    = sub$Name,
    plot      = TRUE,
    plot.opts = list(
      plot.observed   = TRUE,
      plot.model      = TRUE,
      plot.dispersion = FALSE,
      legend          = FALSE,
      lwd             = 1,
      main            = ""
    )
  )
  axis(1, lwd = 2.5, lwd.ticks = 1.2)
  axis(2, lwd = 2.5, lwd.ticks = 1.2)
  box(lwd = 2.5)
  dev.off()
  
  s_proj <- as.data.frame(summary(Nonpareil.set(sub$File, plot = FALSE))) %>%
    mutate(Name = sub$Name, File = sub$File, Project = p) %>%
    select(Project, Name, File, kappa, C, LR, LRstar, modelR, diversity)
  write_csv(s_proj, file.path("nonpareil_per_project", paste0("summary_", p, ".csv")))
}

