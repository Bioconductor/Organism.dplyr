## create views: idmap, genename, pmid(?), ensenbltrans, protein, view_go 
## view_go_all, map, alias, ipi, path

org <- org.Hs.eg.db
conn = dbconn(org.Hs.eg.db)

## use case: map from 'A' to 'B'; expect 1:1 mappings ids = tbl %>% map_ids()
##
## flat table mapping between these identifiers
"ENTREZID",
"ENSEMBL",
"SYMBOL"
"UCSCKG",
"UNIGENE",

sql <- "CREATE TEMPORARY VIEW idmap AS 
        SELECT DISTINCT
            genes._id AS entrez,
            ensembl.ensembl_id AS ensembl, 
            gene_info.symbol AS symbol,
            unigene.unigene_id AS unigene,
            ucsc.ucsc_id AS ucsckkg
        FROM genes 
        JOIN gene_info ON genes._id = gene_info._id
        JOIN ensembl ON genes._id = ensembl._id
        LEFT OUTER JOIN unigene ON genes._id = unigene._id
        LEFT OUTER JOIN ucsc ON genes._id = ucsc._id
        "
dbSendQuery(conn, sql)
idmap <- src_sql("sqlite", conn, path = dbfile(org)) %>% tbl("idmap")

ids = idmap %>% filter(symbol %in% c("BRCA1", "BRCA2")) %>% 
    select(entrez, symbol, unigene, map) %>% distinct

## from ENSEMBL, I want...; ids %>% info(GENENAME)
##
## SQL mapping from central id (ENSEMBL) to...

"ENTREZID",
"GENENAME",
"OMIM",
"ACCNUM",
"REFSEQ",

sql <- "CREATE TEMPORARY VIEW genename AS 
        SELECT DISTINCT
            genes._id AS entrez, 
            gene_info.gene_name AS genename, 
            omim.omim_id AS omim,
            accessions.accession AS accnum, 
            refseq.accession AS refseq,
            pubmed.pubmed_id AS pmid
        FROM genes 
        JOIN gene_info ON genes._id = gene_info._id
        JOIN omim ON genes._id = omim._id
        LEFT OUTER JOIN accessions ON genes._id = accessions._id
        LEFT OUTER JOIN refseq ON genes._id = refseq._id
            AND refseq.accession = accessions.accession
        LEFT OUTER JOIN pubmed ON genes._id = pubmed._id
        "

dbSendQuery(conn, sql)
genename <- src_sql("sqlite", conn, path = dbfile(org)) %>% tbl("genename")

genename %>% filter(entrez %in% c(554, 556)) %>% 
    select(entrez, omim, accnum) %>% distinct

genename %>% filter(entrez %in% c(554, 556)) %>% 
    select(entrez, omim) %>% distinct


## from ENSEMBL, I want...; ids %>% pubmed
##
## SQL mapping from central id (ENSEMBL) to...

"ENTREZID",
"PMID",

sql <- "CREATE TEMPORARY VIEW pmid AS 
        SELECT DISTINCT
            genes._id AS entrez, 
            pubmed.pubmed_id AS pmid
        FROM genes 
        LEFT OUTER JOIN pubmed ON genes._id = pubmed._id
        "

dbSendQuery(conn, sql)
pmid <- src_sql("sqlite", conn, path = dbfile(org)) %>% tbl("pmid")

pmid %>% filter(entrez %in% c(554, 556)) %>% 
    select(entrez, pmid) %>% distinct

## from identifer to transcript: ids %>% transcript_ids()
"ENTREZID",
"ENSEMBLTRANS",

sql <- "CREATE TEMPORARY VIEW ensenbltrans AS 
        SELECT DISTINCT
            genes._id AS entrez, 
            ensembl_trans.trans_id AS ensenbltrans
        FROM genes 
        LEFT OUTER JOIN ensembl_trans ON genes._id = ensembl_trans._id
        "

dbSendQuery(conn, sql)
ensenbltrans <- src_sql("sqlite", conn, path = dbfile(org)) %>% tbl("ensenbltrans")

ensenbltrans %>% filter(entrez %in% c(554, 556)) %>% 
    select(entrez, ensenbltrans) %>% distinct

## from identifer to protien: ids %>% protein_ids()
"ENTREZID",
"ENSEMBLPROT",
"UNIPROT"
"ENZYME",

sql <- "CREATE TEMPORARY VIEW protein AS
        SELECT DISTINCT 
            genes._id AS entrez, 
            ec.ec_number AS enzyme,
            ensembl_prot.prot_id AS ensemblprot,
            uniprot.uniprot_id AS uniprot
        FROM genes 
        LEFT OUTER JOIN ec ON genes._id = ec._id
        LEFT OUTER JOIN uniprot ON genes._id = uniprot._id
        LEFT OUTER JOIN ensembl_prot ON genes._id = ensembl_prot._id
        "
dbSendQuery(conn, sql)
protein <- src_sql("sqlite", conn, path = dbfile(org)) %>% tbl("protein")

protein %>% filter(entrez %in% c(554, 556)) %>% 
    select(entrez, enzyme, uniprot) %>% distinct


## from identifer to protien: ids %>% go_ids()
## mygo = go %>% filter(ontology=="BP") %>% select(id, goid, symbol) %>% distinct
##
"ENTREZID",
"GO",
"EVIDENCE",
"ONTOLOGY",
"GOALL",
"EVIDENCEALL",
"ONTOLOGYALL",

sql <- "CREATE TEMPORARY VIEW view_go AS
        SELECT DISTINCT 
            genes._id AS entrez, 
            go.go_id AS go,
            go.evidence AS evidence, 
            go.ontology AS ontology
        FROM genes 
        LEFT OUTER JOIN go ON genes._id = go._id
        "
dbSendQuery(conn, sql)
view_go <- src_sql("sqlite", conn, path = dbfile(org)) %>% tbl("view_go")

view_go %>% filter(entrez %in% c(554, 556)) %>% 
    select(entrez, go, evidence, ontology) %>% distinct

sql <- "CREATE TEMPORARY VIEW view_go_all AS
        SELECT DISTINCT 
            genes._id AS entrez, 
            go_all.go_id AS goall,
            go_all.evidence AS evidenceall, 
            go_all.ontology AS ontologyall
        FROM genes 
        LEFT OUTER JOIN go_all ON genes._id = go_all._id
        "
dbSendQuery(conn, sql)
view_go_all <- src_sql("sqlite", conn, path = dbfile(org)) %>% tbl("view_go_all")

view_go_all %>% filter(entrez %in% c(554, 556)) %>% 
    select(entrez, goall, evidenceall, ontologyall) %>% distinct



## from identifer to chromosome bands: ids %>% map()
"ENTREZID",
"MAP",

sql <- "CREATE TEMPORARY VIEW map AS 
         SELECT DISTINCT
            genes._id AS entrez, 
            chromosomes.chromosome AS chromosome, 
            chromosome_locations.start_location, 
            chromosome_locations.end_location, 
            chrlengths.length, 
            cytogenetic_location AS map
        FROM genes 
        JOIN chromosomes ON genes._id = chromosomes._id
        JOIN chromosome_locations ON genes._id = chromosome_locations._id 
        JOIN chrlengths ON chromosomes.chromosome = chrlengths.chromosome
        JOIN cytogenetic_locations ON genes._id = cytogenetic_locations._id
        "
dbSendQuery(conn, sql)
map <- src_sql("sqlite", conn, path = dbfile(org)) %>% tbl("map")

map %>% filter(entrez %in% c(554, 556)) %>% 
    select(entrez, chromosome, map) %>% distinct


## "ALIAS"
"ENTREZID",
"ALIAS"
sql <- "CREATE TEMPORARY VIEW alias AS 
         SELECT DISTINCT
            genes._id AS entrez, 
            alias.alias_symbol AS alias
        FROM genes 
        LEFT OUTER JOIN alias ON genes._id = alias._id
        "

##
"ENTREZID",
"IPI"
"PFAM"
"PROSITE"

sql <- "CREATE TEMPORARY VIEW ipi AS
        SELECT DISTINCT 
            genes._id AS entrez, 
            pfam.ipi_id AS ipi,
            pfam.pfam_id AS pfam, 
            prosite.prosite_id AS prosite
        FROM genes 
        LEFT OUTER JOIN pfam ON genes._id = pfam._id
        LEFT OUTER JOIN prosite ON genes._id = prosite._id
            AND pfam.ipi_id = prosite.ipi_id
        "
dbSendQuery(conn, sql)
ipi <- src_sql("sqlite", conn, path = dbfile(org)) %>% tbl("ipi")

ipi %>% filter(entrez %in% c(554, 556)) %>% 
    select(entrez, ipi, pfam, prosite) %>% distinct


##
"ENTREZID",
"PATH"
sql <- "CREATE TEMPORARY VIEW path AS 
        SELECT DISTINCT
            genes._id AS entrez, 
            kegg.path_id AS path
        FROM genes 
        LEFT OUTER JOIN kegg ON genes._id = kegg._id
        "

dbSendQuery(conn, sql)
path <- src_sql("sqlite", conn, path = dbfile(org)) %>% tbl("path")

path %>% filter(entrez %in% c(554, 556)) %>% 
    select(entrez, path) %>% distinct


