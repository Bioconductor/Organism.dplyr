CREATE TABLE IF NOT EXISTS ranges_gene AS
SELECT DISTINCT
    tx_chrom AS gene_chrom,
    MIN(tx_start) AS gene_start,
    MAX(tx_end) AS gene_end,
    tx_strand AS gene_strand,
    gene.gene_id AS entrez
FROM txdb_entrez.gene
JOIN txdb_entrez.transcript ON transcript._tx_id = gene._tx_id
GROUP BY entrez, tx_chrom, tx_strand;

CREATE INDEX IF NOT EXISTS entrez_ranges_gene on ranges_gene (entrez); 

CREATE TABLE IF NOT EXISTS ranges_tx AS
SELECT DISTINCT
    tx_chrom,
    tx_strand,
    tx_start,
    tx_end,
    gene.gene_id AS entrez,
    transcript._tx_id AS tx_id,
    tx_name,
    tx_type
FROM txdb_entrez.transcript
LEFT OUTER JOIN txdb_entrez.gene ON transcript._tx_id = gene._tx_id;

CREATE INDEX IF NOT EXISTS entrez_ranges_tx on ranges_tx (entrez); 

CREATE INDEX IF NOT EXISTS txid_ranges_tx on ranges_tx (tx_id);

CREATE TABLE IF NOT EXISTS ranges_exon AS
SELECT DISTINCT
    exon_chrom,
    exon_strand,
    exon_start,
    exon_end,
    gene.gene_id AS entrez,
    splicing._tx_id AS tx_id,
    exon._exon_id AS exon_id,
    exon_name,
    exon_rank
FROM txdb_entrez.exon
JOIN txdb_entrez.splicing ON exon._exon_id = splicing._exon_id
JOIN txdb_entrez.transcript ON splicing._tx_id = transcript._tx_id
LEFT OUTER JOIN txdb_entrez.gene ON transcript._tx_id = gene._tx_id;

CREATE INDEX IF NOT EXISTS entrez_ranges_exon on ranges_exon (entrez); 

CREATE INDEX IF NOT EXISTS txid_ranges_exon on ranges_exon (tx_id);

CREATE INDEX IF NOT EXISTS exonid_ranges_exon on ranges_exon (exon_id);

CREATE TABLE IF NOT EXISTS ranges_cds AS
SELECT DISTINCT
    cds_chrom,
    cds_strand,
    cds_start,
    cds_end,
    gene.gene_id AS entrez,
    splicing._tx_id AS tx_id,
    cds._cds_id AS cds_id,
    cds_name,
    exon_rank
FROM txdb_entrez.cds
JOIN txdb_entrez.splicing ON cds._cds_id = splicing._cds_id
JOIN txdb_entrez.transcript ON splicing._tx_id = transcript._tx_id
LEFT OUTER JOIN txdb_entrez.gene ON transcript._tx_id = gene._tx_id;

CREATE INDEX IF NOT EXISTS entrez_ranges_cds on ranges_cds (entrez); 

CREATE INDEX IF NOT EXISTS txid_ranges_cds on ranges_cds (tx_id);

CREATE INDEX IF NOT EXISTS cdsid_ranges_cds on ranges_cds (cds_id);

CREATE TABLE IF NOT EXISTS metadata_txdb AS
SELECT * FROM metadata;
