test_that("GBIF concurrency-limit messages are recognised", {
  expect_true(
    .bf_is_gbif_download_limit_error(
      simpleError("A download limitation is exceeded: user has too many simultaneous downloads")
    )
  )
  expect_false(
    .bf_is_gbif_download_limit_error(
      simpleError("Invalid credentials")
    )
  )
})

test_that("account-level slot guard waits until an active job clears", {
  poll_i <- 0L
  sleep_i <- 0L

  fake_list <- function(user, pwd, limit = 1000) {
    poll_i <<- poll_i + 1L
    if (poll_i <= 2L) {
      list(results = data.frame(status = "RUNNING"))
    } else {
      list(results = data.frame(status = "SUCCEEDED"))
    }
  }

  fake_sleep <- function(seconds) {
    sleep_i <<- sleep_i + 1L
    invisible(NULL)
  }

  expect_true(
    isTRUE(
      .bf_wait_for_gbif_download_slot(
        user = "user",
        pwd = "pwd",
        max_active_downloads = 1L,
        poll_seconds = 0,
        max_polls = 5L,
        quiet = TRUE,
        list_fun = fake_list,
        sleep_fun = fake_sleep
      )
    )
  )

  expect_equal(poll_i, 3L)
  expect_equal(sleep_i, 2L)
})

test_that("a GBIF slot-limit rejection retries the same submission", {
  submit_i <- 0L

  fake_list <- function(user, pwd, limit = 1000) {
    list(results = data.frame(status = character()))
  }

  fake_submit <- function(...) {
    submit_i <<- submit_i + 1L

    if (submit_i <= 2L) {
      stop(
        "A download limitation is exceeded: user has too many simultaneous downloads",
        call. = FALSE
      )
    }

    list(key = "000TEST-GBIF-KEY")
  }

  fake_sleep <- function(seconds) invisible(NULL)

  res <- .bf_submit_gbif_with_retry(
    predicate_args = list("fake_predicate"),
    user = "user",
    pwd = "pwd",
    email = "user@example.com",
    label = "test request",
    max_active_downloads = 1L,
    slot_poll_seconds = 0,
    slot_max_polls = 1L,
    submission_max_tries = 5L,
    quiet = TRUE,
    submit_fun = fake_submit,
    list_fun = fake_list,
    sleep_fun = fake_sleep
  )

  expect_equal(res$status, "submitted")
  expect_equal(res$key, "000TEST-GBIF-KEY")
  expect_equal(res$attempts, 3L)
  expect_equal(submit_i, 3L)
})
