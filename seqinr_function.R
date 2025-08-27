# cargo paquete de seqinr
library(seqinr)

# creao una funcion con un solo argumento. Se le debe pasar un archivo fasta como argumento
length_seq <- function(fasta_file) {
  # lee archivo multifasta
  seqs <- read.fasta(fasta_file)
  
  # Crear tabla con resultados
  results <- data.frame(
    Sequence_ID = names(seqs),
    Length = sapply(seqs, length),
    stringsAsFactors = FALSE
  )
  
  return(results)
}

