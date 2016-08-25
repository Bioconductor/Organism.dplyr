CREATE TEMPORARY VIEW IF NOT EXISTS ensenbltrans AS
SELECT DISTINCT
    genes.gene_id AS entrez,
    ensembl_trans.trans_id AS ensenbltrans
FROM genes
LEFT OUTER JOIN ensembl_trans ON genes._id = ensembl_trans._id
