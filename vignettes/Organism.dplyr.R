## ---- echo=FALSE---------------------------------------------------------
suppressPackageStartupMessages({
    library(Organism.dplyr)
    library(GenomicRanges)
    library(ggplot2)
})

## ---- eval=FALSE---------------------------------------------------------
## library(Organism.dplyr)

## ---- eval=FALSE---------------------------------------------------------
## src <- src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene")

## ---- eval=FALSE---------------------------------------------------------
## path <- "path/to/my.sqlite"
## src <- src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene", path)

## ------------------------------------------------------------------------
supportedOrganisms()

## ---- eval=FALSE---------------------------------------------------------
## src <- src_ucsc("human", path)

## ------------------------------------------------------------------------
src <- src_organism(dbpath = hg38light())
src

## ------------------------------------------------------------------------
src_tbls(src)

## ------------------------------------------------------------------------
tbl(src, "id")

## ------------------------------------------------------------------------
colnames(tbl(src, "id"))

## ------------------------------------------------------------------------
tbl(src, "id") %>%
    filter(symbol %like% "SNORD%") %>%
    dplyr::select(entrez, map, ensembl, symbol) %>%
    distinct() %>% arrange(symbol) %>% collect()

## ------------------------------------------------------------------------
inner_join(tbl(src, "id"), tbl(src, "id_go")) %>%
    filter(symbol == "ADA") %>%
    dplyr::select(entrez, ensembl, symbol, go, evidence, ontology)

## ------------------------------------------------------------------------
txcount <- inner_join(tbl(src, "id"), tbl(src, "ranges_tx")) %>%
    dplyr::select(symbol, tx_id) %>%
    group_by(symbol) %>%
    summarise(count = count(symbol)) %>%
    dplyr::select(symbol, count) %>%
    arrange(desc(count)) %>%
    collect(n=Inf)

txcount

library(ggplot2)
ggplot(txcount, aes(x = symbol, y = count)) + 
    geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Transcript count") +
    labs(x = "Symbol") +
    labs(y = "Count")

## ------------------------------------------------------------------------
inner_join(tbl(src, "id"), tbl(src, "ranges_gene")) %>%
    filter(symbol %in% c("ADA", "NAT2")) %>%
    dplyr::select(gene_chrom, gene_start, gene_end, gene_strand,
                  symbol, map) %>%
    collect() %>% GenomicRanges::GRanges()

## ------------------------------------------------------------------------
keytypes(src)

## ------------------------------------------------------------------------
columns(src)

## ------------------------------------------------------------------------
head(keys(src))

## ------------------------------------------------------------------------
head(keys(src, "symbol"))

## ------------------------------------------------------------------------
keytype <- "symbol"
keys <- c("ADA", "NAT2")
columns <- c("entrez", "tx_id", "tx_name","exon_id")
select_tbl(src, keys, columns, keytype)

## ------------------------------------------------------------------------
mapIds(src, keys, column = "tx_name", keytype)
mapIds(src, keys, column = "tx_name", keytype, multiVals="CharacterList")

## ------------------------------------------------------------------------
supportedFilters()

## ------------------------------------------------------------------------
EnsemblFilter("ENSG00000196839")
SymbolFilter("SNORD", "startsWith")

## ------------------------------------------------------------------------
filters <- list(SymbolFilter("SNORD", "startsWith"))
transcripts_tbl(src, filter=filters)
filters <- list(
    SymbolFilter("SNORD", "startsWith"),
    GRangesFilter(GenomicRanges::GRanges("chr15:25062333-25065121"))
)
transcripts(src, filter=filters)

## ------------------------------------------------------------------------
transcripts_tbl(src, filter = list(
    SymbolFilter("ADA"),
    TxStartFilter(44619810,"<")
))

## ------------------------------------------------------------------------
sessionInfo()

