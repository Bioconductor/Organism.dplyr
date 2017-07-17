library(microbenchmark)

fib1 <- function(seq) {
	if (seq == 1) return(1)
    if (seq == 2) return(2)
    return(fib1(seq - 1) + fib1(seq - 2))
}

fib2 <- local({
    memo <- c(1, 1, rep(NA, 100))
    f <- function(x) {
        if(x == 0) return(0)
        if(x < 0) return(NA)
        if(x > length(memo))
        stop("’x’ too big for implementation")
        if(!is.na(memo[x])) return(memo[x])
        ans <- f(x-2) + f(x-1)
        memo[x] <<- ans
        ans
    }
})

x <- sample(20, 100, replace=T)
