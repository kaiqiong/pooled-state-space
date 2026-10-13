## Final design without optimization, for a preprocessing check and simulations.
w <- readRDS("meal_windows_final_65min.rds")
bio <- read.csv(file.path(APP_DATA, "bio.csv"), check.names = FALSE)
bio$subject <- sprintf("CGMacros-%03d", bio$subject)
bio$group <- cut(bio[["A1c PDL (Lab)"]], breaks = c(-Inf, 5.7, 6.5, Inf),
                 labels = c("healthy", "prediabetic", "diabetic"), right = FALSE)
s <- bio[match(w$meta$subject, bio$subject), ]
stopifnot(!anyNA(s[, c("Age", "Gender", "BMI", "group")]))
Xmat <- cbind(Intercept = 1, Age = as.numeric(scale(s$Age)), Sex = as.numeric(s$Gender == "F"),
               BMI = as.numeric(scale(s$BMI)), Prediabetic = as.numeric(s$group == "prediabetic"),
               Diabetic = as.numeric(s$group == "diabetic"))
design <- list(Y = w$Y - w$meta$baseline_gl, Xmat = Xmat,
               meta = cbind(w$meta, s[, c("Age", "Gender", "BMI", "group")]),
               flags = w$flags, slot_time = w$slot_time)
saveRDS(design, "application_design.rds")
cohort <- data.frame(duration_minutes = c(65, 120), windows = c(nrow(w$Y), NA),
                     participants = c(length(unique(w$meta$subject)), NA), slots = c(ncol(w$Y), NA))
w120 <- readRDS("meal_windows_final_120min.rds")
cohort[2, 2:4] <- c(nrow(w120$Y), length(unique(w120$meta$subject)), ncol(w120$Y))
write.csv(cohort, "cohort_summary.csv", row.names = FALSE)
write.csv(data.frame(duration_minutes = rep(c(65, 120), each = length(w$flow)),
                     step = rep(names(w$flow), 2), count = c(w$flow, w120$flow)),
          "cohort_flow.csv", row.names = FALSE)
print(cohort)
stopifnot(nrow(w$Y) == 1174L, ncol(w$Y) == 14L, length(unique(w$meta$subject)) == 44L,
          nrow(w120$Y) == 1059L, ncol(w120$Y) == 25L,
          all(is.finite(design$Y)), all(is.finite(Xmat)))
cat("Final manuscript cohort check passed. Design saved without running fits.\n")
