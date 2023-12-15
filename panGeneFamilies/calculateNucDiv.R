library(pegas)

args <- commandArgs(trailingOnly = TRUE)



readr::write_lines(nuc.div(read.dna(args[1],format = "fasta")),
            file = args[2])
