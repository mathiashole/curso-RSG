# cargo paquete de seqinr
library(seqinr)

# creao una funcion con un solo argumento. Se le debe pasar un archivo fasta como argumento
length_seq <- function(fasta_file) {
  # lee archivo multifasta
  seqs <- read.fasta(fasta_file, as.string = TRUE, forceDNAtolower = FALSE)
  
  # Crear tabla con resultados
  results <- data.frame(
    Sequence_ID = names(seqs),
    Length = sapply(seqs, length),
    GC_Content = sapply(seqs, function(x) {
      round(GC(x) * 100, 2)  # Porcentaje de GC con 2 decimales
    }),
    stringsAsFactors = FALSE
  )
  
  return(results)
}

