# Portable October 13 replication. Run from the repository root; see docs/SIMULATIONS.txt.
source("scripts/simulation/helpers/clustered_summary.R")
run_summary("posvar", getwd(), "results/simulation/tables", "final")
run_summary("boundary", getwd(), "results/simulation/tables", "final")
# Reformat the regenerated same-200 positive-variance summaries for paired tables.
d <- read.csv("results/simulation/tables/part_q_summary_final_paired.csv")
p <- read.csv("results/simulation/tables/part_q_summary_final_paired_differences.csv")
d$empirical_sd_same200 <- d$mean_se / d$se_empirical_sd_ratio
d$bootstrap_minus_coverage <- NA_real_; d$paired_mcse <- NA_real_
for (i in seq_len(nrow(d))) if (d$method[i] != "bootstrap") {
 z <- p[p$tau == d$tau[i] & p$parameter == d$parameter[i] & p$comparison == paste0("bootstrap_minus_", d$method[i]), ]
 stopifnot(nrow(z) == 1L)
 d$bootstrap_minus_coverage[i] <- z$coverage_difference; d$paired_mcse[i] <- z$paired_mcse
}
write.csv(d, "results/simulation/tables/paired_posvar200_comparison.csv", row.names = FALSE)
