CREATE TEMPORARY TABLE IF NOT EXISTS ranges_gene AS
SELECT DISTINCT
    gene.gene_id AS entrez,
    tx_chrom AS chrom,
    gene._tx_id AS txid,
    MIN(tx_start) AS start,
    MAX(tx_end) AS end,
    tx_strand AS strand
FROM txdb_entrez.gene
LEFT OUTER JOIN txdb_entrez.transcript ON txdb_entrez.transcript._tx_id = txdb_entrez.gene._tx_id
GROUP BY entrez, chrom, txid, strand;

CREATE INDEX entrez_ranges_gene on ranges_gene (entrez); 
CREATE INDEX txid_ranges_gene on ranges_gene (txid); 

CREATE TEMPORARY TABLE IF NOT EXISTS ranges_tx AS
SELECT DISTINCT
    gene.gene_id AS entrez,
    transcript._tx_id AS txid,
    transcript.tx_chrom AS chrom,
    transcript.tx_strand AS strand,
    transcript.tx_start AS start,
    transcript.tx_end AS end
FROM txdb_entrez.gene
LEFT OUTER JOIN txdb_entrez.transcript ON txdb_entrez.transcript._tx_id = txdb_entrez.gene._tx_id;

CREATE INDEX entrez_ranges_tx on ranges_tx (entrez); 
CREATE INDEX txid_ranges_tx on ranges_tx (txid);

CREATE TEMPORARY TABLE IF NOT EXISTS ranges_exon AS
SELECT DISTINCT
    gene.gene_id AS entrez,
    transcript._tx_id AS txid,
    exon._exon_id AS exonid,
    exon.exon_chrom AS chrom,
    exon.exon_strand AS strand,
    exon.exon_start AS start,
    exon.exon_end AS end,
    splicing.exon_rank AS rank
FROM txdb_entrez.gene
LEFT OUTER JOIN txdb_entrez.transcript ON txdb_entrez.transcript._tx_id = txdb_entrez.gene._tx_id
LEFT OUTER JOIN txdb_entrez.splicing ON txdb_entrez.splicing._tx_id = txdb_entrez.gene._tx_id
LEFT OUTER JOIN txdb_entrez.exon ON txdb_entrez.exon._exon_id = txdb_entrez.splicing._exon_id;

CREATE INDEX entrez_ranges_exon on ranges_exon (entrez); 
CREATE INDEX txid_ranges_exon on ranges_exon (txid);
CREATE INDEX exonid_ranges_exon on ranges_exon (exonid);

CREATE TEMPORARY TABLE IF NOT EXISTS ranges_cds AS
SELECT DISTINCT
    gene.gene_id AS entrez,
    transcript._tx_id AS txid,
    cds._cds_id AS cdsid,
    cds.cds_chrom AS chrom,
    cds.cds_strand AS strand,
    cds.cds_start AS start,
    cds.cds_end AS end
FROM txdb_entrez.gene
LEFT OUTER JOIN txdb_entrez.transcript ON txdb_entrez.transcript._tx_id = txdb_entrez.gene._tx_id
LEFT OUTER JOIN txdb_entrez.splicing ON txdb_entrez.splicing._tx_id = txdb_entrez.gene._tx_id
LEFT OUTER JOIN txdb_entrez.cds ON txdb_entrez.cds._cds_id = txdb_entrez.splicing._cds_id;

CREATE INDEX entrez_ranges_cds on ranges_cds (entrez); 
CREATE INDEX txid_ranges_cds on ranges_cds (txid);
CREATE INDEX cdsid_ranges_cds on ranges_cds (cdsid);

CREATE TEMPORARY TABLE IF NOT EXISTS id_ranges AS
SELECT DISTINCT
    gene.gene_id AS entrez,
    transcript._tx_id AS txid,
    transcript.tx_name AS txname,
    transcript.tx_type AS txtype,
    exon._exon_id AS exonid,
    exon.exon_name AS exonname,
    cds._cds_id AS cdsid,
    cds.cds_name AS cdsname
FROM txdb_entrez.gene
LEFT OUTER JOIN txdb_entrez.transcript ON txdb_entrez.transcript._tx_id = txdb_entrez.gene._tx_id
LEFT OUTER JOIN txdb_entrez.splicing ON txdb_entrez.splicing._tx_id = txdb_entrez.gene._tx_id
LEFT OUTER JOIN txdb_entrez.exon ON txdb_entrez.exon._exon_id = txdb_entrez.splicing._exon_id
LEFT OUTER JOIN txdb_entrez.cds ON txdb_entrez.cds._cds_id = txdb_entrez.splicing._cds_id;

CREATE INDEX entrez_id_ranges on id_ranges (entrez); 
CREATE INDEX txid_id_ranges on id_ranges (txid); 
CREATE INDEX exonid_id_ranges on id_ranges (exonid);
CREATE INDEX cdsid_id_ranges on id_ranges (cdsid);
