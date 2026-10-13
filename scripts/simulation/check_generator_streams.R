# Bounded reconstruction check: only function definitions are evaluated.
# Each dataset and B=999 resampling stream has its own seed. Therefore omitting
# the superseded B99/B199 calculations cannot change any final stream.
get_generator <- function(path, environment) {
  expressions <- as.list(parse(path))
  candidates <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                         identical(x[[2]], as.name("gen")), expressions)
  stopifnot(length(candidates) == 1L)
  eval(candidates[[1]], environment)
  environment$gen
}
env <- new.env(parent = globalenv())
env$Xmat <- cbind(Intercept = 1, Age = rep(c(-1, 0, 0.5, 1), each = 2))
env$subj <- rep(c("p1", "p2", "p3", "p4"), each = 2)
env$N <- nrow(env$Xmat); env$Tt <- 5L; env$G <- 4L
with(env, {
  beta0 <- c(0.4, 0.2); a1_0 <- 2; P1_0 <- 1.5
  se2_0 <- 0.8; r0 <- 400; sp2_0 <- 0
})
boundary <- get_generator("scripts/simulation/clustered_boundary.R", env)
fixed <- get_generator("scripts/simulation/clustered_fixed_variance.R", env)
free <- get_generator("scripts/simulation/clustered_free_variance.R", env)
checks <- list()
record <- function(name, pass) checks[[length(checks) + 1L]] <<- data.frame(check = name, pass = pass)
for (tau in c(0, 0.4)) for (m in c(1L, 200L, 500L, 1000L)) {
  seed <- 410000 + 10000 * round(10 * tau) + m
  set.seed(seed); q <- fixed(tau); end_q <- .Random.seed
  set.seed(seed); r <- free(tau); end_r <- .Random.seed
  record(sprintf("fixed_free_same_dataset_tau%s_m%d", tau, m), identical(q, r) && identical(end_q, end_r))
  for (which in c("boundary", "positive")) {
    generator <- if (which == "boundary") boundary else fixed
    base <- if (which == "boundary") 900000 else 410000
    ds <- base + 10000 * round(10 * tau) + m
    set.seed(ds); y <- generator(tau)
    invisible(rnorm(123)) # arbitrary preceding pilot/bootstrap stream
    set.seed(ds); y_again <- generator(tau)
    record(sprintf("dataset_reseed_%s_tau%s_m%d", which, tau, m), identical(y, y_again))
    bs <- ds + 5000000
    draw <- function() lapply(seq_len(999L), function(b) sample(unique(env$subj), env$G, TRUE))
    set.seed(bs); ids <- draw()
    invisible(rnorm(257)); set.seed(bs); ids_again <- draw()
    record(sprintf("bootstrap_reseed_%s_tau%s_m%d", which, tau, m), identical(ids, ids_again))
  }
}
out <- do.call(rbind, checks)
dir.create("results/simulation/validation", recursive = TRUE, showWarnings = FALSE)
write.csv(out, "results/simulation/validation/generator_stream_checks.csv", row.names = FALSE)
stopifnot(all(out$pass))
cat(nrow(out), "generator and bootstrap-stream checks passed\n")
