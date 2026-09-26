testthat::test_that(
  "central downloader quotes Windows system-curl path and URL arguments",
  {
    download_fun <- getFromNamespace(
      "bf_download_cached",
      "biofetchR"
    )
    
    txt <- paste(
      deparse(
        body(download_fun),
        width.cutoff = 500L
      ),
      collapse = "\n"
    )
    
    # Command-line curl output paths must be quoted.
    testthat::expect_match(
      txt,
      "shQuote(tmp)",
      fixed = TRUE
    )
    
    # URLs passed to command-line curl must also be quoted.
    testthat::expect_match(
      txt,
      "shQuote(url)",
      fixed = TRUE
    )
    
    # The old unsafe forms must no longer exist.
    testthat::expect_false(
      grepl(
        '"-o", tmp',
        txt,
        fixed = TRUE
      )
    )
    
    testthat::expect_false(
      grepl(
        "args <- c(args, url)",
        txt,
        fixed = TRUE
      )
    )
  }
)

testthat::test_that(
  "central downloader accepts a destination path containing spaces and a standalone hyphen",
  {
    
    root <- file.path(
      tempdir(),
      "OneDrive - Queen's University Belfast",
      "biofetchR curl path test"
    )
    
    dir.create(
      root,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    # Keep the synthetic source URL simple. The regression target is the
    # destination/cache path containing spaces and a standalone hyphen.
    source_file <- tempfile(
      pattern = "biofetchR_source_",
      fileext = ".bin"
    )
    
    destination <- file.path(
      root,
      "downloaded file.bin"
    )
    
    payload <- as.raw(
      rep(
        c(65L, 66L, 67L, 68L),
        2048L
      )
    )
    
    writeBin(
      payload,
      source_file
    )
    
    source_norm <- normalizePath(
      source_file,
      winslash = "/",
      mustWork = TRUE
    )
    
    file_url <- if (.Platform$OS.type == "windows") {
      paste0(
        "file:///",
        source_norm
      )
    } else {
      paste0(
        "file://",
        source_norm
      )
    }
    
    out <- bf_download_cached(
      url = file_url,
      dest = destination,
      force_refresh = TRUE,
      quiet = TRUE,
      retries = 0L,
      min_bytes = 1024L,
      validate_not_html = FALSE,
      resume = TRUE,
      keep_partial = TRUE
    )
    
    testthat::expect_true(
      file.exists(out)
    )
    
    testthat::expect_equal(
      unname(file.info(out)$size),
      unname(file.info(source_file)$size)
    )
    
    downloaded <- readBin(
      out,
      what = "raw",
      n = file.info(out)$size
    )
    
    testthat::expect_identical(
      downloaded,
      payload
    )
  }
)