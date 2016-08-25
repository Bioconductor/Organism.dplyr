CREATE TEMPORARY VIEW IF NOT EXISTS idmap AS 
SELECT DISTINCT
    genes.gene_id AS entrez,
    cytogenetic_location AS map, 
    ensembl.ensembl_id AS ensembl, 
    gene_info.symbol AS symbol,
    unigene.unigene_id AS unigene
FROM genes 
JOIN gene_info ON genes._id = gene_info._id
JOIN cytogenetic_locations ON genes._id = cytogenetic_locations._id
JOIN ensembl ON genes._id = ensembl._id
LEFT OUTER JOIN unigene ON genes._id = unigene._id
