# Replace with your own chosen phenotype and image ID
phenotype <- "Actin dots"
image_id <- "19465"
# Read the Fiji Results CSV (no R packages needed)
d <- read.csv("fiji_data.csv", check.names = FALSE)
needed <- c("Area", "Mean", "Circ.", "AR")
if (!all(needed %in% names(d))) {
  stop("Missing Fiji columns. Recheck Set Measurements.")
}
d <- d[complete.cases(d[, needed]), ]
if (nrow(d) < 5) stop("Measure five clear cells.")
d$Cell <- paste0("Cell ", seq_len(nrow(d)))
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
# Save numerical summary
s <- data.frame(phenotype = phenotype, image_id = image_id,
                cells_measured = nrow(d), mean_area_px2 = mean(d$Area),
                mean_circularity = mean(d[["Circ."]]),
                mean_aspect_ratio = mean(d$AR),
                mean_actin_signal = mean(d$Mean))
write.csv(s, "results/cell_summary.csv", row.names = FALSE)
print(s)
# Figure 1: Cell area versus circularity
png("figures/area_vs_circularity.png", width = 900, height = 650)
plot(d$Area, d[["Circ."]], pch = 19, col = "steelblue",
     xlab = "Cell area (pixels squared)", ylab = "Circularity (0 to 1)",
     main = paste(phenotype, "-", image_id))
text(d$Area, d[["Circ."]], labels = d$Cell, pos = 3, cex = 0.75)
dev.off()
# Figure 2: Relative actin signal of measured cells
png("figures/actin_intensity.png", width = 900, height = 650)
barplot(d$Mean, names.arg = d$Cell, col = "steelblue",
        ylab = "Mean actin fluorescence (arbitrary units)",
        main = paste("Actin signal -", image_id))
dev.off()

# FIGURE 2: MYH6 exon-intron variant schematic

# Exon 13 ends at transcript coordinate 1477.
# The variant affects the first intronic nucleotide
# immediately following the exon.

exons <- data.frame(
  exon = c("Exon 13", "Exon 14"),
  start = c(0, 70),
  end = c(30, 100)
)

fig2 <- ggplot() +
  
  # Exon boxes
  geom_rect(
    data = exons,
    aes(
      xmin = start,
      xmax = end,
      ymin = 0.4,
      ymax = 0.7
    ),
    fill = "steelblue",
    color = "black"
  ) +
  
  # Intron connecting the exons
  geom_segment(
    aes(
      x = 30,
      xend = 70,
      y = 0.55,
      yend = 0.55
    ),
    linewidth = 1.5,
    color = "grey50"
  ) +
  
  # Splice-donor variant marker
  geom_point(
    aes(x = 31, y = 0.55),
    color = "red",
    size = 4
  ) +
  
  # Exon labels
  annotate(
    "text",
    x = 15,
    y = 0.85,
    label = "Exon 13",
    size = 5
  ) +
  
  annotate(
    "text",
    x = 85,
    y = 0.85,
    label = "Exon 14",
    size = 5
  ) +
  
  # Intron label
  annotate(
    "text",
    x = 50,
    y = 0.32,
    label = "Intron (not to scale)",
    size = 4
  ) +
  
  # Variant label
  annotate(
    "label",
    x = 43,
    y = 1.2,
    label = "c.1410+1G>A\nSplice-donor variant",
    color = "red",
    size = 4
  ) +
  
  labs(
    title = "MYH6 Splice-Site Variant",
    subtitle = paste(
      "Reference transcript:", transcript_accession,
      "| ClinVar ID:", clinvar_id
    ),
    x = "Schematic representation of exon-intron structure",
    y = NULL
  ) +
  
  coord_cartesian(
    xlim = c(-5, 105),
    ylim = c(0.2, 1.5)
  ) +
  
  theme_bw() +
  
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 16)
  )

print(fig2)

ggsave(
  "figures/02_variant_position.png",
  fig2,
  width = 10,
  height = 5,
  dpi = 300
)
