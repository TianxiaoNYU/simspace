# install_r_deps.R

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}
renv::restore()