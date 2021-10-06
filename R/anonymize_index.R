#clean up the csv files by deleting the first 3 columns

library("purrr")

dir.create(file.path("data","Index_anon"))
files <- list.files(file.path("data","Index"), pattern="\\.csv$", full.names=TRUE)
walk(files, function(x) {
  read_csv(x)[,-c(1:3)] %>% 
  write_csv(file.path("data","Index_anon", basename(x)))
})
