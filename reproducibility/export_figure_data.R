# Export the numerical data underlying Figures 1-3 and 5-7.
invisible(Sys.setlocale('LC_ALL', 'English_United States.utf8'))
args <- commandArgs(trailingOnly = FALSE)
root <- dirname(normalizePath(sub('^--file=', '', grep('^--file=', args, value = TRUE)[1]), winslash = '/'))
d <- read.csv(file.path(root, 'final_analysis_data.csv'), check.names = FALSE)
stopifnot(nrow(d) == 1003, !anyDuplicated(d$case_id))
out <- file.path(dirname(root), 'figure_data')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
save_csv <- function(x, name) write.csv(x, file.path(out, name), row.names = FALSE, na = '', fileEncoding = 'UTF-8')
distribution <- function(variable, codes, labels) {
  n <- sum(!is.na(d[[variable]]))
  count <- vapply(codes, function(k) sum(d[[variable]] == k, na.rm = TRUE), numeric(1))
  data.frame(code = codes, label = labels, count = count, denominator = n, percent = 100 * count / n)
}
save_csv(distribution('time', c(2, 3, 1), c('06:00-18:00', '18:00-24:00', '00:00-06:00')), 'Figure_1_time.csv')
save_csv(distribution('road', c(1, 5, 3, 2, 4), c('Urban roads', 'Rural roads', 'National/provincial roads', 'Expressways', 'Residential/parking')), 'Figure_1_road.csv')
age <- data.frame(age_code = 1:4, age_group = c('24 years or younger', '25 to 44 years', '45 to 64 years', '65 years or older'))
for (role in 1:2) for (sex in 1:2) {
  field <- paste(c('initiating', 'affected')[role], c('male', 'female')[sex], sep = '_')
  age[[field]] <- vapply(1:4, function(k) sum(d[[paste0('age', role)]] == k & d[[paste0('sex', role)]] == sex, na.rm = TRUE), numeric(1))
}
age$origin_initiating_male <- -age$initiating_male
age$origin_initiating_female <- -age$initiating_female
stopifnot(sum(age$initiating_male + age$initiating_female) == 794,
          sum(age$affected_male + age$affected_female) == 487, age$initiating_female[4] == 0)
save_csv(age, 'Figure_2_joint_age_sex.csv')
for (variable in c('cause', 'process')) {
  codes <- seq_len(if (variable == 'cause') 7 else 6)
  result <- do.call(rbind, lapply(1:2, function(s) {
    observed <- !is.na(d$sex1) & d$sex1 == s & !is.na(d[[variable]])
    n <- sum(observed)
    count <- vapply(codes, function(k) sum(observed & d[[variable]] == k, na.rm = TRUE), numeric(1))
    data.frame(sex_code = s, sex = c('Male', 'Female')[s], category = codes, count = count, denominator = n, percent = 100 * count / n)
  }))
  save_csv(result, paste0('Figure_3_', variable, '.csv'))
}
for (variable in c('time', 'road')) {
  codes <- seq_len(if (variable == 'time') 3 else 5)
  labels <- if (variable == 'time') c('00:00-06:00', '06:00-18:00', '18:00-24:00') else c('Urban roads', 'Expressways', 'National/provincial roads', 'Residential/parking', 'Rural roads')
  tb <- table(factor(d[[variable]], levels = codes), factor(d$conflict, levels = 1:5))
  test <- suppressWarnings(chisq.test(tb, correct = FALSE))
  adjusted <- test$stdres
  holm <- matrix(p.adjust(2 * pnorm(-abs(as.vector(adjusted))), method = 'holm'), nrow = nrow(tb))
  result <- do.call(rbind, lapply(seq_along(codes), function(i) {
    data.frame(row_code = codes[i], row_label = labels[i], type = paste0('T', 1:5), count = as.numeric(tb[i, ]),
               adjusted_residual = as.numeric(adjusted[i, ]), pearson_residual = as.numeric(test$residuals[i, ]),
               holm_p = as.numeric(holm[i, ]), significant = as.numeric(holm[i, ]) < .05)
  }))
  save_csv(result, paste0('Figure_', if (variable == 'time') 5 else 6, '_residuals.csv'))
}
weighted <- read.csv(file.path(root, 'tables', 'lca_type_weighted.csv'))
stopifnot(identical(names(weighted), c('type', 'C1', 'C2', 'C3')), all(abs(rowSums(weighted[, -1]) - 1) < 1e-10))
weighted[, -1] <- weighted[, -1] * 100
save_csv(weighted, 'Figure_7_posterior_weighted_percent.csv')
cat('Exported 8 figure-data CSVs.\n')
