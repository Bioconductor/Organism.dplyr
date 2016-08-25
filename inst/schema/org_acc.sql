CREATE TEMPORARY VIEW IF NOT EXISTS acc AS 
SELECT DISTINCT
    genes.gene_id AS entrez, 
    accessions.accession AS accnum, 
    refseq.accession AS refseq
FROM genes 
JOIN gene_info ON genes._id = gene_info._id
LEFT OUTER JOIN accessions ON genes._id = accessions._id
LEFT OUTER JOIN refseq ON genes._id = refseq._id
    AND refseq.accession = accessions.accession
