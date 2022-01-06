fp1 <- system.file("extdata", "trend_gads_2020.db", package = "eatGADS")
fp2 <- system.file("extdata", "trend_gads_2015.db", package = "eatGADS")
fp3 <- system.file("extdata", "trend_gads_2010.db", package = "eatGADS")
load(system.file("extdata", "linking_error.rda", package = "eatRep"))
#leDF <- data.frame()

test_that("3 mp trend", {
  expect_message(checkLEs(c(fp1, fp2, fp3), lErr),
               "The following variables have linking errors but are not variables in data base 1: 'value', 'valueTransfBista'")
  
  #le_mes <- capture_messages(out <- checkLEs(filePaths = c(fp1, fp2, fp3), leDF = leDF))
  #expect_equal(out$dep_notIn_nam[[3]], "transfBista")
  #expect_equal(out$key_notIn_nam[[3]], "parameter")
})
