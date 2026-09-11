# Calculate probe weighting factor
# Translation of TransQuant matlab script

library(Biostrings)
library(tidyverse)
# library(scattermore)

compute_transquant <- function(target_fasta, probe_delim){
  transcript_fasta <- readDNAStringSet(target_fasta, format = "fasta")
  transcript_sequence <- as.character(transcript_fasta)
  L <- str_length(transcript_sequence)
  
  str_revcomp <- function(seq){
    chartr("ATGC", "TACG", seq) %>% stringi::stri_reverse()
  }
  
  probes <- read_delim(probe_delim, delim = "\t", col_names = FALSE) %>%
    setNames("probe_sequence") %>%
    mutate(probe_sequence = str_to_upper(probe_sequence)) %>%
    mutate(probe_sequence = str_revcomp(probe_sequence)) %>%
    mutate(probe_length = str_length(probe_sequence)) %>%
    mutate(
      start = str_locate(transcript_sequence, probe_sequence)[,1],
      end = str_locate(transcript_sequence, probe_sequence)[,2]
    ) %>%
    mutate(mid_point = (start + end) / 2)
  
  N <- list()
  for (i in 1:L){
    probe_count_till_here <- probes %>%
      filter(mid_point < i) %>%
      nrow()
    N[[i]] <- probe_count_till_here
  }
  
  N_df <- map_dfr(N, ~ tibble(count = .x), .id = "coordinate") %>%
    mutate(proportion = (count / nrow(probes)) * 100) %>%
    mutate(coordinate = as.integer(coordinate))
  
  W <- (1 / L) * sum(N_df$count) / N[[L]]
  
  plot <- N_df %>%
    ggplot(aes(x = coordinate, y = proportion)) +
    geom_point(size = 0.1, alpha = 0.5) +
    labs(title = paste0("Transcript length = ", L, "bp & W = ", W),
         x = "Coordinate along the transcript", 
         y = "Proportion of bound probes") + 
    theme_minimal()
  
  output <- list(
    "W" = W,
    "L" = L,
    "N" = N,
    "plot" = plot
  )
  
  return(output)
}

# x <- compute_transquant("~/Desktop/Seq.txt", "~/Desktop/Probes.txt")

# ## Target transcript
# transcript_fasta <- readDNAStringSet("~/Desktop/Seq.txt", format = "fasta")
# transcript_sequence <- as.character(transcript_fasta)
# L <- str_length(transcript_sequence)
# 
# ## smFISH probes (location is 1-based)
# str_revcomp <- function(seq){
#   chartr("ATGC", "TACG", seq) %>% stringi::stri_reverse()
# }
# 
# probes <- read_delim("~/Desktop/Probes.txt", delim = "\t", col_names = FALSE) %>%
#   setNames("probe_sequence") %>%
#   mutate(probe_sequence = str_to_upper(probe_sequence)) %>%
#   mutate(probe_sequence = str_revcomp(probe_sequence)) %>%
#   mutate(probe_length = str_length(probe_sequence)) %>%
#   mutate(
#     start = str_locate(transcript_sequence, probe_sequence)[,1],
#     end = str_locate(transcript_sequence, probe_sequence)[,2]
#   ) %>%
#   mutate(mid_point = (start + end) / 2)
# 
# ## Calculate N(x) which us the number of probes found until coordinate x
# N <- list()
# for (i in 1:L){
#   probe_count_till_here <- probes %>%
#     filter(mid_point < i) %>%
#     nrow()
#   N[[i]] <- probe_count_till_here
# }
# 
# N_df <- map_dfr(N, ~ tibble(count = .x), .id = "coordinate") %>%
#   mutate(proportion = (count / nrow(probes)) * 100) %>%
#   mutate(coordinate = as.integer(coordinate))
# 
# ## Calculate W
# W <- (1 / L) * sum(N_df$count) / N[[L]]
# 
# ## Plot probe coverage across transcript
# plot <- N_df %>%
#   ggplot(aes(x = coordinate, y = proportion)) +
#   geom_point(size = 0.1, alpha = 0.5) +
#   # geom_scattermore() + 
#   labs(title = paste0("Transcript length = ", L, "bp & W = ", W),
#        x = "Coordinate along the transcript", 
#        y = "Proportion of bound probes") + 
#   theme_minimal()
# 
# output <- list(
#   "W" = W,
#   "L" = L,
#   "N" = N,
#   "plot" = plot
# )
# 
# output$plot