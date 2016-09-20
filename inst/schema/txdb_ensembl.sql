CREATE TEMPORARY VIEW IF NOT EXISTS ranges_gene AS
SELECT DISTINCT
    gene.gene_id AS ensembl,
    tx_chrom AS chrom,
    MIN(tx_start) AS start,
    MAX(tx_end) AS end,
    tx_strand AS strand
FROM txdb_ensembl.gene
LEFT OUTER JOIN txdb_ensembl.transcript ON txdb_ensembl.transcript._tx_id = txdb_ensembl.gene._tx_id
GROUP BY txdb_ensembl.gene.gene_id;

CREATE TEMPORARY VIEW IF NOT EXISTS ranges_tx AS
SELECT DISTINCT
    gene.gene_id AS ensembl,
    transcript._tx_id AS txid,
    transcript.tx_chrom AS chrom,
    transcript.tx_strand AS strand,
    transcript.tx_start AS start,
    transcript.tx_end AS end
FROM txdb_ensembl.gene
LEFT OUTER JOIN txdb_ensembl.transcript ON txdb_ensembl.transcript._tx_id = txdb_ensembl.gene._tx_id;

CREATE TEMPORARY VIEW IF NOT EXISTS ranges_exon AS
SELECT DISTINCT
    gene.gene_id AS ensembl,
    transcript._tx_id AS txid,
    exon._exon_id AS exonid,
    exon.exon_chrom AS chrom,
    exon.exon_strand AS strand,
    exon.exon_start AS start,
    exon.exon_end AS end,
    splicing.exon_rank AS rank
FROM txdb_ensembl.gene
LEFT OUTER JOIN txdb_ensembl.transcript ON txdb_ensembl.transcript._tx_id = txdb_ensembl.gene._tx_id
LEFT OUTER JOIN txdb_ensembl.splicing ON txdb_ensembl.splicing._tx_id = txdb_ensembl.gene._tx_id
LEFT OUTER JOIN txdb_ensembl.exon ON txdb_ensembl.exon._exon_id = txdb_ensembl.splicing._exon_id;

CREATE TEMPORARY VIEW IF NOT EXISTS ranges_cds AS
SELECT DISTINCT
    gene.gene_id AS ensembl,
    transcript._tx_id AS txid,
    cds._cds_id AS cdsid,
    cds.cds_chrom AS chrom,
    cds.cds_strand AS strand,
    cds.cds_start AS start,
    cds.cds_end AS end
FROM txdb_ensembl.gene
LEFT OUTER JOIN txdb_ensembl.transcript ON txdb_ensembl.transcript._tx_id = txdb_ensembl.gene._tx_id
LEFT OUTER JOIN txdb_ensembl.splicing ON txdb_ensembl.splicing._tx_id = txdb_ensembl.gene._tx_id
LEFT OUTER JOIN txdb_ensembl.cds ON txdb_ensembl.cds._cds_id = txdb_ensembl.splicing._cds_id;

CREATE TEMPORARY VIEW IF NOT EXISTS id_ranges AS
SELECT DISTINCT
    gene.gene_id AS ensembl,
    transcript._tx_id AS txid,
    transcript.tx_name AS txname,
    transcript.tx_type AS txtype,
    exon._exon_id AS exonid,
    exon.exon_name AS exonname,
    cds._cds_id AS cdsid,
    cds.cds_name AS cdsname
FROM txdb_ensembl.gene
LEFT OUTER JOIN txdb_ensembl.transcript ON txdb_ensembl.transcript._tx_id = txdb_ensembl.gene._tx_id
LEFT OUTER JOIN txdb_ensembl.splicing ON txdb_ensembl.splicing._tx_id = txdb_ensembl.gene._tx_id
LEFT OUTER JOIN txdb_ensembl.exon ON txdb_ensembl.exon._exon_id = txdb_ensembl.splicing._exon_id
LEFT OUTER JOIN txdb_ensembl.cds ON txdb_ensembl.cds._cds_id = txdb_ensembl.splicing._cds_id;
