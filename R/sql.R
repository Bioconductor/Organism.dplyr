##
## create views for org.Hs.eg.db
##
## flat1
sql <- "CREATE TEMPORARY VIEW flat AS 
        SELECT DISTINCT
            genes._id AS id, 
            gene_info.gene_name AS genename, 
            gene_info.symbol AS symbol, 
            chromosomes.chromosome AS chromosome, 
            chromosome_locations.start_location, 
            chromosome_locations.end_location, 
            chrlengths.length, 
            cytogenetic_location AS location, 
            ensembl.ensembl_id AS ensembl, 
            ec.ec_number AS ec, 
            unigene.unigene_id AS unigene
        FROM genes 
        JOIN gene_info ON genes._id = gene_info._id
        JOIN chromosomes ON genes._id = chromosomes._id
        JOIN chromosome_locations ON genes._id = chromosome_locations._id 
        JOIN chrlengths ON chromosomes.chromosome = chrlengths.chromosome
        JOIN cytogenetic_locations ON genes._id = cytogenetic_locations._id
        JOIN ensembl ON genes._id = ensembl._id
        LEFT OUTER JOIN ec ON genes._id = ec._id
        LEFT OUTER JOIN unigene ON genes._id = unigene._id
        "
dbSendQuery(conn, sql)

## flat2
sql <- "CREATE TEMPORARY VIEW flat AS 
        SELECT DISTINCT
            genes._id AS id, 
            gene_info.gene_name AS genename, 
            gene_info.symbol AS symbol, 
            chromosomes.chromosome AS chromosome, 
            cytogenetic_location AS location, 
            ensembl.ensembl_id AS ensembl, 
            ec.ec_number AS ec, 
            unigene.unigene_id AS unigene, 
            ensembl_prot.prot_id AS prot, 
            ensembl_trans.trans_id AS trans
        FROM genes 
        JOIN gene_info ON genes._id = gene_info._id
        JOIN chromosomes ON genes._id = chromosomes._id
        JOIN cytogenetic_locations ON genes._id = cytogenetic_locations._id
        JOIN ensembl ON genes._id = ensembl._id
        LEFT OUTER JOIN ec ON genes._id = ec._id
        LEFT OUTER JOIN unigene ON genes._id = unigene._id
        LEFT OUTER JOIN ensembl_prot ON genes._id = ensembl_prot._id
        LEFT OUTER JOIN ensembl_trans ON genes._id = ensembl_trans._id
        "

## go 
sql <- "CREATE TEMPORARY VIEW view_go AS
        SELECT DISTINCT 
            genes._id AS id, 
            gene_info.symbol AS symbol,
            go.go_id AS goid,
            go.evidence AS evidence, 
            go.ontology AS ontology
        FROM genes 
        JOIN gene_info ON genes._id = gene_info._id
        LEFT OUTER JOIN go ON genes._id = go._id
        "
dbSendQuery(conn, sql)

## protein 
sql <- "CREATE TEMPORARY VIEW protein AS
        SELECT DISTINCT 
            genes._id AS id, 
            gene_info.symbol AS symbol,
            pfam.ipi_id AS ipiid,
            pfam.pfam_id AS pfam, 
            prosite.prosite_id AS prosite
        FROM genes 
        JOIN gene_info ON genes._id = gene_info._id
        LEFT OUTER JOIN pfam ON genes._id = pfam._id
        LEFT OUTER JOIN prosite ON genes._id = prosite._id
        AND pfam.ipi_id = prosite.ipi_id
        "
dbSendQuery(conn, sql)

## dplyr query
org <- org.Hs.eg.db
conn = dbconn(org)
flat <- src_sql("sqlite", conn, path = dbfile(org)) %>% tbl("flat")
go <- src_sql("sqlite", conn, path = dbfile(org)) %>% tbl("view_go")
protein <- src_sql("sqlite", conn, path = dbfile(org)) %>% tbl("protein")

result = flat %>% 
            filter(symbol == "PTEN") %>% 
            select(ensembl, symbol, location, ec, unigene) 
distinct(result)

result = inner_join(flat, go) %>% 
            filter(unigene == "Hs.500466") %>% 
            select(symbol, goid, evidence, ontology)
distinct(result)

result = inner_join(flat, protein) %>% 
            filter(ec == "3.1.3.67") %>% 
            select(symbol, pfam, prosite)
distinct(result)
