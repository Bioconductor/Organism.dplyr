CREATE TEMPORARY VIEW IF NOT EXISTS ensenbltrans AS
SELECT DISTINCT
    genes.gene_id AS entrez,
    unigene.unigene_id AS unigene,
    ensembl_trans.trans_id AS ensenbltrans
FROM genes
LEFT OUTER JOIN unigene ON genes._id = unigene._id
LEFT OUTER JOIN ensembl_trans ON genes._id = ensembl_trans._id
