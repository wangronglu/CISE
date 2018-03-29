context("Input Validation")

test_that("threshold tol should be positive",{
  data(A)
  expect_error(MGRAF1(A=A,K=1,tol=-0.01))
  expect_error(MGRAF2(A=A,K=1,tol=-0.01))
  expect_error(MGRAF3(A=A,K=1,tol=-0.01))
})
