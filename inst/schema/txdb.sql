CREATE TEMPORARY VIEW IF NOT EXISTS ranges_gene AS
SELECT DISTINCT
    gene.gene_id AS entrez,
    tx_chrom AS chrom,
    MIN(tx_start) AS start,
    MAX(tx_end) AS end,
    tx_strand AS strand
FROM txdb.gene
LEFT OUTER JOIN txdb.transcript ON txdb.transcript._tx_id = txdb.gene._tx_id
GROUP BY txdb.gene.gene_id;

CREATE TEMPORARY VIEW IF NOT EXISTS ranges_tx AS
SELECT DISTINCT
    gene.gene_id AS entrez,
    transcript._tx_id AS id,
    transcript.tx_chrom AS chrom,
    transcript.tx_strand AS strand,
    transcript.tx_start AS start,
    transcript.tx_end AS end
FROM txdb.gene
LEFT OUTER JOIN txdb.transcript ON txdb.transcript._tx_id = txdb.gene._tx_id;

CREATE TEMPORARY VIEW IF NOT EXISTS ranges_exon AS
SELECT DISTINCT
    gene.gene_id AS entrez,
    transcript._tx_id AS txid,
    exon._exon_id AS exonid,
    exon.exon_chrom AS exonchrom,
    exon.exon_strand AS exonstrand,
    exon.exon_start AS exonstart,
    exon.exon_end AS exonend,
    splicing.exon_rank AS exonrank
FROM txdb.gene
LEFT OUTER JOIN txdb.transcript ON txdb.transcript._tx_id = txdb.gene._tx_id
LEFT OUTER JOIN txdb.splicing ON txdb.splicing._tx_id = txdb.gene._tx_id
LEFT OUTER JOIN txdb.exon ON txdb.exon._exon_id = txdb.splicing._exon_id;

CREATE TEMPORARY VIEW IF NOT EXISTS ranges_cds AS
SELECT DISTINCT
    gene.gene_id AS entrez,
    transcript._tx_id AS txid,
    cds._cds_id AS cdsid,
    cds.cds_chrom AS cdschrom,
    cds.cds_strand AS cdsstrand,
    cds.cds_start AS cdsstart,
    cds.cds_end AS cdsend
FROM txdb.gene
LEFT OUTER JOIN txdb.transcript ON txdb.transcript._tx_id = txdb.gene._tx_id
LEFT OUTER JOIN txdb.splicing ON txdb.splicing._tx_id = txdb.gene._tx_id
LEFT OUTER JOIN txdb.cds ON txdb.cds._cds_id = txdb.splicing._cds_id;

CREATE TEMPORARY VIEW IF NOT EXISTS id_ranges AS
SELECT DISTINCT
    gene.gene_id AS entrez,
    transcript._tx_id AS txid,
    transcript.tx_name AS txname,
    transcript.tx_type AS txtype,
    exon._exon_id AS exonid,
    exon.exon_name AS exonname,
    cds._cds_id AS cdsid,
    cds.cds_name AS cdsname
FROM txdb.gene
LEFT OUTER JOIN txdb.transcript ON txdb.transcript._tx_id = txdb.gene._tx_id
LEFT OUTER JOIN txdb.splicing ON txdb.splicing._tx_id = txdb.gene._tx_id
LEFT OUTER JOIN txdb.exon ON txdb.exon._exon_id = txdb.splicing._exon_id
LEFT OUTER JOIN txdb.cds ON txdb.cds._cds_id = txdb.splicing._cds_id;
