coale_model <- function(sample_size, recombination = FALSE){
  coale_list <- list(n = sample_size, coale_event=list())
  if (recombination) {
    coale_list <- append(coale_list, list(recomb=list(recomb_rate=recombination, recomb_event=list())))
  }
  coale_list
}
# construct class for the model