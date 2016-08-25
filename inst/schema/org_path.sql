CREATE TEMPORARY VIEW IF NOT EXISTS path AS
SELECT DISTINCT
    genes.gene_id AS entrez,
    kegg.path_id AS path
FROM genes
LEFT OUTER JOIN kegg ON genes._id = kegg._id
