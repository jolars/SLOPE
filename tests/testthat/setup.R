# Capture plots without plotting
dont_plot <- function(x, ...) {
  tmp <- tempfile()
  grDevices::png(tmp)
  p <- plot(x, ...)
  grDevices::dev.off()
  unlink(tmp)
  invisible(p)
}

# Suppress printning
dont_print <- function(x, ...) {
  utils::capture.output(y <- print(x, ...))
  invisible(y)
}
