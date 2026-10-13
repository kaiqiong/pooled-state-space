args <- commandArgs(trailingOnly=TRUE)
required <- c('numDeriv','statmod','nlme','ggplot2','patchwork','dplyr','tidyr')
missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if (length(missing) && '--install' %in% args) {
  install.packages(missing,repos='https://cloud.r-project.org')
  missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
}
cat(R.version.string,'\n')
for (p in setdiff(required,missing)) cat(p,as.character(packageVersion(p)),'\n')
if (length(missing)) stop('Missing R packages: ',paste(missing,collapse=', '),'. Run Rscript scripts/check_environment.R --install')
cat('R dependencies available. Tested versions are recorded in config/R-packages-tested.csv.\n')
