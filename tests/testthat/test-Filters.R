context("Filters")

test_that("BasicFilter", {
    expect_error(SymbolFilter())
    expect_error(SymbolFilter(1, ">"))
    expect_error(SymbolFilter(1, "foo"))
    expect_error(Tx_endFilter("foo"))
    expect_error(Tx_endFilter(c(1,2), ">"))
    expect_error(SymbolFilter(c("foo","bar"), "startsWith"))
})

test_that("GRangesFilter", {
    expect_error(GRangesFilter())
    expect_error(GRangesFilter("foo"))
})
