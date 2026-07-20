# Route all graphics from plot_it = TRUE tests to a null device so that no
# Rplots.pdf (or any file) is produced while running the test suite.
grDevices::pdf(NULL)
