# Route all graphics from the tests to a null device so that no Rplots.pdf (or
# any file) is produced while running the test suite. The device is closed
# once the suite ends.
withr::local_pdf(NULL, .local_envir = testthat::teardown_env())
