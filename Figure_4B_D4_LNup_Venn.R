# Day 4 LN-upregulated genes shared across Basalt, Oakley and F₁
# Venn diagram based on the seven overlap counts
# No input file required
# Outputs: PDF, 600 dpi PNG and CSV with region counts

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos = "https://cloud.r-project.org")
}
library(ggplot2)

# Counts for the seven mutually exclusive regions
regions <- data.frame(
  region = c("Basalt only", "Oakley only", "F1 only", "Basalt & Oakley only",
             "Basalt & F1 only", "Oakley & F1 only", "All three"),
  count = c(98L, 6272L, 12871L, 21L, 30L, 1590L, 1L),
  stringsAsFactors = FALSE
)

with(as.list(setNames(regions$count, c("b", "o", "f", "bo", "bf", "of", "all"))), {
  stopifnot(b + bo + bf + all == 150L,
            o + bo + of + all == 7884L,
            f + bf + of + all == 14492L,
            of + all == 1591L)
})

# Circles are the same size; their areas do not represent gene counts
circle <- function(cx, cy, radius, name, n = 900L) {
  theta <- seq(0, 2*pi, length.out = n)
  data.frame(x = cx + radius*cos(theta), y = cy + radius*sin(theta), group = name)
}

r <- 1.65
circles <- rbind(circle(-0.90, 0.55, r, "Basalt"),
                 circle( 0.90, 0.55, r, "Oakley"),
                 circle( 0.00,-0.80, r, "F1"))
circles$group <- factor(circles$group, levels = c("Basalt", "Oakley", "F1"))

numbers <- data.frame(
  x = c(-1.65, 1.65, 0, 0, -0.94, 0.94, 0),
  y = c( 0.85, 0.85,-1.85, 1.30,-0.60,-0.60,0.03),
  label = format(regions$count, big.mark = ",", scientific = FALSE, trim = TRUE)
)

# Convert font size from points to millimetres for ggplot2
pt_to_mm <- function(pt) pt * 25.4 / 72.27
font_pt <- 14

p <- ggplot() +
  geom_polygon(data = circles,
               aes(x = x, y = y, group = group, fill = group),
               alpha = 0.52, colour = "#28323A", linewidth = 0.55) +
  geom_text(data = numbers, aes(x = x, y = y, label = label),
            size = pt_to_mm(font_pt), colour = "#151515", family = "sans") +
  annotate("text", x = -2.20, y = 2.43, label = "Basalt",
           fontface = "bold", size = pt_to_mm(font_pt), family = "sans") +
  annotate("text", x = 2.20, y = 2.43, label = "Oakley",
           fontface = "bold", size = pt_to_mm(font_pt), family = "sans") +
  annotate("text", x = 0, y = -2.65, label = "F\u2081",
           fontface = "bold", size = pt_to_mm(font_pt), family = "sans") +
  scale_fill_manual(values = c(Basalt = "#D9D9D9",
                               Oakley = "#87B5A3", F1 = "#D7B5D8")) +
  coord_fixed(xlim = c(-2.85, 2.85), ylim = c(-2.92, 2.88),
              expand = FALSE, clip = "off") +
  theme_void(base_size = font_pt) +
  theme(legend.position = "none", plot.margin = margin(12, 12, 12, 12))

# Save beside the open script in RStudio, or in the working directory
script_path <- tryCatch(rstudioapi::getSourceEditorContext()$path,
                        error = function(e) "")
out_dir <- if (length(script_path) == 1L && nzchar(script_path) && file.exists(script_path)) {
  dirname(normalizePath(script_path, winslash = "/"))
} else getwd()

base <- file.path(out_dir, "Figure_3B_D4_LNup_Venn_14PT")
write.csv(regions, paste0(base, "_region_counts.csv"), row.names = FALSE)

ggsave(paste0(base, ".pdf"), plot = p,
       width = 7.2, height = 7.3, units = "in", device = "pdf",
       bg = "white", useDingbats = FALSE)
ggsave(paste0(base, "_600dpi.png"), plot = p,
       width = 7.2, height = 7.3, units = "in", dpi = 600, bg = "white")

print(p)
message("Saved Venn diagram to: ", out_dir,
        "\nOakley and F₁ shared LN-up genes: 1,590 + 1 = 1,591")
