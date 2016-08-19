CREATE TEMPORARY VIEW idmap AS 
SELECT DISTINCT
    genes.gene_id AS entrez,
    ensembl.ensembl_id AS ensembl, 
    gene_info.symbol AS symbol,
    unigene.unigene_id AS unigene,
    ucsc.ucsc_id AS ucsckkg
FROM genes 
JOIN gene_info ON genes._id = gene_info._id
JOIN ensembl ON genes._id = ensembl._id
LEFT OUTER JOIN unigene ON genes._id = unigene._id
LEFT OUTER JOIN ucsc ON genes._id = ucsc._id
