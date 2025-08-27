# cargo paquete de seqinr
library(seqinr)

# creao una funcion con un solo argumento. Se le debe pasar un archivo fasta como argumento
length_seq <- function(fasta_file) {
  # lee archivo multifasta
  seqs <- read.fasta(fasta_file)
  
  # creo una tabla id de cada secuencia y su largo, guardandolo en dos columnas id y length
  results <- data.frame(
    id = names(seqs),
    length = sapply(seqs, length) # 
  )
  
  return(results)
}

