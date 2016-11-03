CREATE TABLE IF NOT EXISTS id_accession AS
SELECT DISTINCT
    genes.gene_id AS entrez,
    accessions.accession AS accnum,
    refseq.accession AS refseq
FROM genes
LEFT OUTER JOIN accessions ON genes._id = accessions._id
LEFT OUTER JOIN refseq ON genes._id = refseq._id
    AND refseq.accession = accessions.accession;

CREATE TABLE IF NOT EXISTS id_transcript AS
SELECT DISTINCT
    genes.gene_id AS entrez,
    unigene.unigene_id AS unigene,
    ensembl_trans.trans_id AS ensembltrans
FROM genes
LEFT OUTER JOIN unigene ON genes._id = unigene._id
LEFT OUTER JOIN ensembl_trans ON genes._id = ensembl_trans._id;

CREATE TABLE IF NOT EXISTS id AS
SELECT DISTINCT
    genes.gene_id AS entrez,
    mgi.mgi_id AS mgi,
    ensembl.ensembl_id AS ensembl,
    gene_info.symbol AS symbol,
    gene_info.gene_name AS genename,
    alias.alias_symbol AS alias
FROM genes
JOIN gene_info ON genes._id = gene_info._id
JOIN mgi ON genes._id = mgi._id
JOIN ensembl ON genes._id = ensembl._id
LEFT OUTER JOIN alias ON genes._id = alias._id;

CREATE INDEX IF NOT EXISTS entrez_id on id (entrez);

CREATE INDEX IF NOT EXISTS ensembl_id on id (ensembl);

CREATE INDEX IF NOT EXISTS symbol_id on id (symbol);

CREATE TABLE IF NOT EXISTS id_pm AS
SELECT DISTINCT
    genes.gene_id AS entrez,
    pubmed.pubmed_id AS pmid
FROM genes
LEFT OUTER JOIN pubmed ON genes._id = pubmed._id;

CREATE TABLE IF NOT EXISTS id_protein AS
SELECT DISTINCT
    genes.gene_id AS entrez,
    ec.ec_number AS enzyme,
    ensembl_prot.prot_id AS ensemblprot,
    uniprot.uniprot_id AS uniprot,
    pfam.ipi_id AS ipi,
    pfam.pfam_id AS pfam,
    prosite.prosite_id AS prosite
FROM genes
LEFT OUTER JOIN ec ON genes._id = ec._id
LEFT OUTER JOIN uniprot ON genes._id = uniprot._id
LEFT OUTER JOIN ensembl_prot ON genes._id = ensembl_prot._id
LEFT OUTER JOIN pfam ON genes._id = pfam._id
LEFT OUTER JOIN prosite ON genes._id = prosite._id
    AND pfam.ipi_id = prosite.ipi_id;

CREATE TABLE IF NOT EXISTS id_go AS
SELECT DISTINCT
    genes.gene_id AS entrez,
    go.go_id AS go,
    go.evidence AS evidence,
    go.ontology AS ontology
FROM genes
LEFT OUTER JOIN go ON genes._id = go._id;

CREATE TABLE IF NOT EXISTS id_go_all AS
SELECT DISTINCT
    genes.gene_id AS entrez,
    go_all.go_id AS goall,
    go_all.evidence AS evidenceall,
    go_all.ontology AS ontologyall
FROM genes
LEFT OUTER JOIN go_all ON genes._id = go_all._id;
