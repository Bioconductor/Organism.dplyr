CREATE TEMPORARY VIEW IF NOT EXISTS txdb AS 
SELECT DISTINCT
    gene.gene_id AS entrez,
    transcript._tx_id AS txid,
    transcript.tx_name AS txname, 
    transcript.tx_type AS txtype, 
    transcript.tx_chrom AS txchrom, 
    transcript.tx_strand AS txstrand, 
    transcript.tx_start AS txstart, 
    transcript.tx_end AS txend, 
    exon._exon_id AS exonid,
    exon.exon_name AS exonname, 
    exon.exon_chrom AS exonchrom, 
    exon.exon_strand AS exonstrand, 
    exon.exon_start AS exonstart, 
    exon.exon_end AS exonend, 
    splicing.exon_rank AS exonrank, 
    cds._cds_id AS cdsid, 
    cds.cds_name AS cdsname, 
    cds.cds_chrom AS cdschrom, 
    cds.cds_strand AS cdsstrand, 
    cds.cds_start AS cdsstart, 
    cds.cds_end AS cdsend
FROM gene
LEFT OUTER JOIN transcript ON transcript._tx_id = gene._tx_id
LEFT OUTER JOIN splicing ON splicing._tx_id = gene._tx_id
LEFT OUTER JOIN exon ON exon._exon_id = splicing._exon_id
LEFT OUTER JOIN cds ON cds._cds_id = splicing._cds_id
