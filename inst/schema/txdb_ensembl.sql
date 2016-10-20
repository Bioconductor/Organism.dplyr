CREATE TEMPORARY TABLE IF NOT EXISTS ranges_gene AS
SELECT DISTINCT
    tx_chrom,
    MIN(tx_start) AS gene_start,
    MAX(tx_end) AS gene_end,
    tx_strand,
    gene.gene_id AS ensembl
FROM txdb_ensembl.gene
LEFT OUTER JOIN txdb_ensembl.transcript ON transcript._tx_id = gene._tx_id
GROUP BY ensembl, tx_chrom, tx_strand;

CREATE INDEX IF NOT EXISTS ensembl_ranges_gene on ranges_gene (ensembl); 

CREATE TEMPORARY TABLE IF NOT EXISTS ranges_tx AS
SELECT DISTINCT
    tx_chrom,
    tx_strand,
    tx_start,
    tx_end,
    gene.gene_id AS ensembl,
    transcript._tx_id AS tx_id,
    tx_name,
    tx_type
FROM txdb_ensembl.transcript
LEFT OUTER JOIN txdb_ensembl.gene ON transcript._tx_id = gene._tx_id;

CREATE INDEX IF NOT EXISTS ensembl_ranges_tx on ranges_tx (ensembl); 

CREATE INDEX IF NOT EXISTS txid_ranges_tx on ranges_tx (tx_id);

CREATE TEMPORARY TABLE IF NOT EXISTS ranges_exon AS
SELECT DISTINCT
    exon_chrom,
    exon_strand,
    exon_start,
    exon_end,
    gene.gene_id AS ensembl,
    splicing._tx_id AS tx_id,
    exon._exon_id AS exon_id,
    exon_name,
    exon_rank
FROM txdb_ensembl.splicing
LEFT OUTER JOIN txdb_ensembl.exon ON exon._exon_id = splicing._exon_id
LEFT OUTER JOIN txdb_ensembl.transcript ON splicing._tx_id = transcript._tx_id
LEFT OUTER JOIN txdb_ensembl.gene ON transcript._tx_id = gene._tx_id;

CREATE INDEX IF NOT EXISTS ensembl_ranges_exon on ranges_exon (ensembl); 

CREATE INDEX IF NOT EXISTS txid_ranges_exon on ranges_exon (tx_id);

CREATE INDEX IF NOT EXISTS exonid_ranges_exon on ranges_exon (exon_id);

CREATE TEMPORARY TABLE IF NOT EXISTS ranges_cds AS
SELECT DISTINCT
    cds_chrom,
    cds_strand,
    cds_start,
    cds_end,
    gene.gene_id AS ensembl,
    splicing._tx_id AS tx_id,
    cds._cds_id AS cds_id,
    cds_name,
    exon_rank
FROM txdb_ensembl.splicing
LEFT OUTER JOIN txdb_ensembl.cds ON cds._cds_id = splicing._cds_id
LEFT OUTER JOIN txdb_ensembl.transcript ON splicing._tx_id = transcript._tx_id
LEFT OUTER JOIN txdb_ensembl.gene ON transcript._tx_id = gene._tx_id
WHERE splicing._tx_id NOT IN 
    (SELECT _tx_id FROM (SELECT distinct _tx_id, _cds_id from splicing)  
         GROUP BY _tx_id HAVING count(*) = 1 AND _cds_id is NULL);

CREATE INDEX IF NOT EXISTS ensembl_ranges_cds on ranges_cds (ensembl); 

CREATE INDEX IF NOT EXISTS txid_ranges_cds on ranges_cds (tx_id);

CREATE INDEX IF NOT EXISTS cdsid_ranges_cds on ranges_cds (cds_id);
