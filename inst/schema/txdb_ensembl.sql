CREATE TEMPORARY TABLE IF NOT EXISTS ranges_gene AS
SELECT DISTINCT
    gene.gene_id AS ensenbl,
    tx_chrom AS chrom,
    gene._tx_id AS txid,
    MIN(tx_start) AS start,
    MAX(tx_end) AS end,
    tx_strand AS strand
FROM txdb_ensenbl.gene
LEFT OUTER JOIN txdb_ensenbl.transcript ON txdb_ensenbl.transcript._tx_id = txdb_ensenbl.gene._tx_id
GROUP BY ensenbl, chrom, txid, strand;

CREATE INDEX ensenbl_ranges_gene on ranges_gene (ensenbl); 

CREATE INDEX txid_ranges_gene on ranges_gene (txid); 

CREATE TEMPORARY TABLE IF NOT EXISTS ranges_tx AS
SELECT DISTINCT
    gene.gene_id AS ensenbl,
    transcript._tx_id AS txid,
    transcript.tx_chrom AS chrom,
    transcript.tx_strand AS strand,
    transcript.tx_start AS start,
    transcript.tx_end AS end,
    transcript.tx_name AS txname,
    transcript.tx_type AS txtype
FROM txdb_ensenbl.gene
LEFT OUTER JOIN txdb_ensenbl.transcript ON txdb_ensenbl.transcript._tx_id = txdb_ensenbl.gene._tx_id;

CREATE INDEX ensenbl_ranges_tx on ranges_tx (ensenbl); 

CREATE INDEX txid_ranges_tx on ranges_tx (txid);

CREATE TEMPORARY TABLE IF NOT EXISTS ranges_exon AS
SELECT DISTINCT
    gene.gene_id AS ensenbl,
    transcript._tx_id AS txid,
    exon._exon_id AS exonid,
    exon.exon_chrom AS chrom,
    exon.exon_strand AS strand,
    exon.exon_start AS start,
    exon.exon_end AS end,
    splicing.exon_rank AS rank,
    exon.exon_name AS exonname
FROM txdb_ensenbl.gene
LEFT OUTER JOIN txdb_ensenbl.transcript ON txdb_ensenbl.transcript._tx_id = txdb_ensenbl.gene._tx_id
LEFT OUTER JOIN txdb_ensenbl.splicing ON txdb_ensenbl.splicing._tx_id = txdb_ensenbl.gene._tx_id
LEFT OUTER JOIN txdb_ensenbl.exon ON txdb_ensenbl.exon._exon_id = txdb_ensenbl.splicing._exon_id;

CREATE INDEX ensenbl_ranges_exon on ranges_exon (ensenbl); 

CREATE INDEX txid_ranges_exon on ranges_exon (txid);

CREATE INDEX exonid_ranges_exon on ranges_exon (exonid);

CREATE TEMPORARY TABLE IF NOT EXISTS ranges_cds AS
SELECT DISTINCT
    gene.gene_id AS ensenbl,
    transcript._tx_id AS txid,
    cds._cds_id AS cdsid,
    cds.cds_chrom AS chrom,
    cds.cds_strand AS strand,
    cds.cds_start AS start,
    cds.cds_end AS end,
    cds.cds_name AS cdsname
FROM txdb_ensenbl.gene
LEFT OUTER JOIN txdb_ensenbl.transcript ON txdb_ensenbl.transcript._tx_id = txdb_ensenbl.gene._tx_id
LEFT OUTER JOIN txdb_ensenbl.splicing ON txdb_ensenbl.splicing._tx_id = txdb_ensenbl.gene._tx_id
LEFT OUTER JOIN txdb_ensenbl.cds ON txdb_ensenbl.cds._cds_id = txdb_ensenbl.splicing._cds_id;

CREATE INDEX ensenbl_ranges_cds on ranges_cds (ensenbl); 

CREATE INDEX txid_ranges_cds on ranges_cds (txid);

CREATE INDEX cdsid_ranges_cds on ranges_cds (cdsid);
