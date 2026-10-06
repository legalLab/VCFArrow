# Save a ggplot as <base>.pdf, <base>.svg and <base>.png (6 x 4 inches,
# transparent background), as the assess_vcf_*() functions do.
#
# 'svglite' is only suggested: SVG files are written with it when it is
# installed, otherwise with R's own Cairo SVG device (grDevices::svg). If
# neither is available the SVG file is skipped with a message.

.save_plot <- function(plt, base) {
  save <- function(ext, device) {
    ggplot2::ggsave(
      plt,
      filename = paste0(base, ".", ext),
      device = device,
      width = 6,
      height = 4,
      bg = "transparent",
      limitsize = FALSE
    )
  }
  save("pdf", "pdf")
  if (requireNamespace("svglite", quietly = TRUE)) {
    save("svg", "svg")
  } else if (isTRUE(capabilities("cairo"))) {
    save("svg", grDevices::svg)
  } else {
    cli::cli_inform(c(
      "i" = "SVG plot not written: install {.pkg svglite} to write SVG files."
    ))
  }
  save("png", "png")
  invisible(NULL)
}
