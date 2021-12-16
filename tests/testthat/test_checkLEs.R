fp1 <- system.file("extdata", "trend_gads_2020.db", package = "eatGADS")
fp2 <- system.file("extdata", "trend_gads_2015.db", package = "eatGADS")
fp3 <- system.file("extdata", "trend_gads_2010.db", package = "eatGADS")
#leDF <- data.frame()

test_that("3 mp trend", {
  expect_error(checkLEs(filePaths = c(fp1, fp2, fp3), leDF = mtcars), 
               "Incorrect variable names in 'leDF'.")
  
  #le_mes <- capture_messages(out <- checkLEs(filePaths = c(fp1, fp2, fp3), leDF = leDF))
  #expect_equal(out$dep_notIn_nam[[3]], "transfBista")
  #expect_equal(out$key_notIn_nam[[3]], "parameter")
})
