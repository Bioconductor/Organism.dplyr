CREATE TEMPORARY VIEW IF NOT EXISTS view_alias AS
SELECT DISTINCT
    genes.gene_id AS entrez,
    alias.alias_symbol AS alias
FROM genes
LEFT OUTER JOIN alias ON genes._id = alias._id
