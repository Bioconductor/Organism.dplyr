CREATE TABLE IF NOT EXISTS id_accession AS
SELECT DISTINCT
    genes.gene_id AS entrez,
    ensembl.ensembl_id AS ensembl,
    accessions.accession AS accnum,
    refseq.accession AS refseq
FROM genes
JOIN accessions ON genes._id = accessions._id
LEFT OUTER JOIN refseq ON genes._id = refseq._id
    AND refseq.accession = accessions.accession
LEFT OUTER JOIN ensembl ON genes._id = ensembl._id;
    
CREATE INDEX IF NOT EXISTS entrez_accession ON id_accession (entrez);

CREATE INDEX IF NOT EXISTS ensembl_accession ON id_accession (ensembl);

CREATE TABLE IF NOT EXISTS id_transcript AS
SELECT DISTINCT
    genes.gene_id AS entrez,
    ensembl.ensembl_id AS ensembl,
    unigene.unigene_id AS unigene,
    ensembl_trans.trans_id AS ensembltrans
FROM genes
LEFT OUTER JOIN ensembl ON genes._id = ensembl._id
LEFT OUTER JOIN unigene ON genes._id = unigene._id
LEFT OUTER JOIN ensembl_trans ON genes._id = ensembl_trans._id;

CREATE INDEX IF NOT EXISTS entrez_transcript ON id_transcript (entrez);

CREATE INDEX IF NOT EXISTS ensembl_transcript ON id_transcript (ensembl);

CREATE TABLE IF NOT EXISTS id AS
SELECT DISTINCT
    genes.gene_id AS entrez,
    ensembl.ensembl_id AS ensembl,
    cytogenetic_location AS map,
    gene_info.symbol AS symbol,
    gene_info.gene_name AS genename,
    alias.alias_symbol AS alias
FROM genes
JOIN gene_info ON genes._id = gene_info._id
LEFT OUTER JOIN ensembl ON genes._id = ensembl._id
LEFT OUTER JOIN cytogenetic_locations ON genes._id = cytogenetic_locations._id
JOIN alias ON genes._id = alias._id;

CREATE INDEX IF NOT EXISTS entrez_id ON id (entrez);

CREATE INDEX IF NOT EXISTS ensembl_id ON id (ensembl);

CREATE INDEX IF NOT EXISTS symbol_id ON id (symbol);

CREATE TABLE IF NOT EXISTS id_flybase AS
SELECT DISTINCT
    genes.gene_id AS entrez,
    ensembl.ensembl_id AS ensembl,
    flybase.flybase_id AS flybase,
    flybase_cg.flybase_cg_id AS flybase_cg,
    flybase_prot.prot_id AS flybase_prot
FROM genes 
LEFT OUTER JOIN ensembl ON genes._id = ensembl._id
LEFT OUTER JOIN flybase ON genes._id = flybase._id
LEFT OUTER JOIN flybase_cg ON genes._id = flybase_cg._id
LEFT OUTER JOIN flybase_prot ON genes._id = flybase_prot._id;

CREATE INDEX IF NOT EXISTS entrez_flybase ON id_flybase (entrez);

CREATE INDEX IF NOT EXISTS ensembl_flybase ON id_flybase (ensembl);

CREATE INDEX IF NOT EXISTS flybase_flybase ON id_flybase (flybase);

CREATE INDEX IF NOT EXISTS flybase_cg_flybase ON id_flybase (flybase_cg);

CREATE TABLE IF NOT EXISTS id_pm AS
SELECT DISTINCT
    genes.gene_id AS entrez,
    ensembl.ensembl_id AS ensembl,
    pubmed.pubmed_id AS pmid
FROM genes
JOIN pubmed ON genes._id = pubmed._id
LEFT OUTER JOIN ensembl ON genes._id = ensembl._id;

CREATE INDEX IF NOT EXISTS entrez_pm ON id_pm (entrez);

CREATE INDEX IF NOT EXISTS ensembl_pm ON id_pm (ensembl);

CREATE TABLE IF NOT EXISTS id_protein AS
SELECT DISTINCT
    genes.gene_id AS entrez,
    ensembl.ensembl_id AS ensembl,
    ec.ec_number AS enzyme,
    ensembl_prot.prot_id AS ensemblprot,
    uniprot.uniprot_id AS uniprot
FROM genes
LEFT OUTER JOIN ensembl ON genes._id = ensembl._id
LEFT OUTER JOIN ec ON genes._id = ec._id
LEFT OUTER JOIN uniprot ON genes._id = uniprot._id
LEFT OUTER JOIN ensembl_prot ON genes._id = ensembl_prot._id;
    
CREATE INDEX IF NOT EXISTS entrez_protein ON id_protein (entrez);

CREATE INDEX IF NOT EXISTS ensembl_protein ON id_protein (ensembl);

CREATE TABLE IF NOT EXISTS id_go AS
SELECT DISTINCT
    genes.gene_id AS entrez,
    ensembl.ensembl_id AS ensembl,
    go.go_id AS go,
    go.evidence AS evidence,
    go.ontology AS ontology
FROM genes
JOIN go ON genes._id = go._id
LEFT OUTER JOIN ensembl ON genes._id = ensembl._id;

CREATE INDEX IF NOT EXISTS entrez_go ON id_go (entrez);

CREATE INDEX IF NOT EXISTS ensembl_go ON id_go (ensembl);

CREATE TABLE IF NOT EXISTS id_go_all AS
SELECT DISTINCT
    genes.gene_id AS entrez,
    ensembl.ensembl_id AS ensembl,
    go_all.go_id AS goall,
    go_all.evidence AS evidenceall,
    go_all.ontology AS ontologyall
FROM genes
JOIN go_all ON genes._id = go_all._id
LEFT OUTER JOIN ensembl ON genes._id = ensembl._id;

CREATE INDEX IF NOT EXISTS entrez_go_all ON id_go_all (entrez);

CREATE INDEX IF NOT EXISTS ensembl_go_all ON id_go_all (ensembl);

CREATE TABLE IF NOT EXISTS metadata_org AS
SELECT * FROM metadata;
