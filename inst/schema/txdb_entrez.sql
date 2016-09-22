CREATE TEMPORARY TABLE IF NOT EXISTS ranges_gene AS
SELECT DISTINCT
    gene.gene_id AS entrez,
    gene._tx_id AS tx_id,
    transcript.tx_chrom AS tx_chrom,
    MIN(tx_start) AS gene_start,
    MAX(tx_end) AS gene_end,
    transcript.tx_strand AS tx_strand,
    transcript.tx_name AS tx_name
FROM txdb_entrez.gene
LEFT OUTER JOIN txdb_entrez.transcript ON txdb_entrez.transcript._tx_id = txdb_entrez.gene._tx_id
GROUP BY entrez, chrom, txid, strand, txname;

CREATE INDEX IF NOT EXISTS entrez_ranges_gene on ranges_gene (entrez); 

CREATE INDEX IF NOT EXISTS txid_ranges_gene on ranges_gene (txid); 

CREATE TEMPORARY TABLE IF NOT EXISTS ranges_tx AS
SELECT DISTINCT
    gene.gene_id AS entrez,
    transcript._tx_id AS txid,
    transcript.tx_chrom AS chrom,
    transcript.tx_strand AS strand,
    transcript.tx_start AS txstart,
    transcript.tx_end AS txend,
    transcript.tx_name AS txname,
    transcript.tx_type AS txtype
FROM txdb_entrez.gene
LEFT OUTER JOIN txdb_entrez.transcript ON txdb_entrez.transcript._tx_id = txdb_entrez.gene._tx_id;

CREATE INDEX IF NOT EXISTS entrez_ranges_tx on ranges_tx (entrez); 

CREATE INDEX IF NOT EXISTS txid_ranges_tx on ranges_tx (txid);

CREATE TEMPORARY TABLE IF NOT EXISTS ranges_exon AS
SELECT DISTINCT
    gene.gene_id AS entrez,
    transcript._tx_id AS txid,
    exon._exon_id AS exonid,
    exon.exon_chrom AS chrom,
    exon.exon_strand AS strand,
    exon.exon_start AS exonstart,
    exon.exon_end AS exonend,
    exon.exon_name AS exonname,
    splicing.exon_rank AS rank
FROM txdb_entrez.gene
LEFT OUTER JOIN txdb_entrez.transcript ON txdb_entrez.transcript._tx_id = txdb_entrez.gene._tx_id
LEFT OUTER JOIN txdb_entrez.splicing ON txdb_entrez.splicing._tx_id = txdb_entrez.gene._tx_id
LEFT OUTER JOIN txdb_entrez.exon ON txdb_entrez.exon._exon_id = txdb_entrez.splicing._exon_id;

CREATE INDEX IF NOT EXISTS entrez_ranges_exon on ranges_exon (entrez); 

CREATE INDEX IF NOT EXISTS txid_ranges_exon on ranges_exon (txid);

CREATE INDEX IF NOT EXISTS exonid_ranges_exon on ranges_exon (exonid);

CREATE TEMPORARY TABLE IF NOT EXISTS ranges_cds AS
SELECT DISTINCT
    gene.gene_id AS entrez,
    transcript._tx_id AS txid,
    cds._cds_id AS cdsid,
    cds.cds_chrom AS chrom,
    cds.cds_strand AS strand,
    cds.cds_start AS cdsstart,
    cds.cds_end AS cdsend,
    cds.cds_name AS cdsname,
    splicing.exon_rank AS rank
FROM txdb_entrez.gene
LEFT OUTER JOIN txdb_entrez.transcript ON txdb_entrez.transcript._tx_id = txdb_entrez.gene._tx_id
LEFT OUTER JOIN txdb_entrez.splicing ON txdb_entrez.splicing._tx_id = txdb_entrez.gene._tx_id
LEFT OUTER JOIN txdb_entrez.cds ON txdb_entrez.cds._cds_id = txdb_entrez.splicing._cds_id;

CREATE INDEX IF NOT EXISTS entrez_ranges_cds on ranges_cds (entrez); 

CREATE INDEX IF NOT EXISTS txid_ranges_cds on ranges_cds (txid);

CREATE INDEX IF NOT EXISTS cdsid_ranges_cds on ranges_cds (cdsid);
