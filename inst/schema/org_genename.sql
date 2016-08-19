CREATE TEMPORARY VIEW genename AS 
SELECT DISTINCT
    genes.gene_id AS entrez, 
    gene_info.gene_name AS genename, 
    omim.omim_id AS omim,
    accessions.accession AS accnum, 
    refseq.accession AS refseq
FROM genes 
JOIN gene_info ON genes._id = gene_info._id
JOIN omim ON genes._id = omim._id
LEFT OUTER JOIN accessions ON genes._id = accessions._id
LEFT OUTER JOIN refseq ON genes._id = refseq._id
    AND refseq.accession = accessions.accession
