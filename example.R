devtools::load_all()

library('TxDb.Hsapiens.UCSC.hg38.knownGene')

src <- src_organism('TxDb.Hsapiens.UCSC.hg38.knownGene')

smbl <- SymbolFilter("ADA")
