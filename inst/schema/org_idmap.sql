CREATE TEMPORARY VIEW IF NOT EXISTS idmap AS 
SELECT DISTINCT
    genes.gene_id AS entrez,
    cytogenetic_location AS map, 
    ensembl.ensembl_id AS ensembl, 
    gene_info.symbol AS symbol,
    gene_info.gene_name AS genename, 
    alias.alias_symbol AS alias
FROM genes 
JOIN gene_info ON genes._id = gene_info._id
JOIN cytogenetic_locations ON genes._id = cytogenetic_locations._id
JOIN ensembl ON genes._id = ensembl._id
LEFT OUTER JOIN alias ON genes._id = alias._id
