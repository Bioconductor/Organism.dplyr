CREATE TEMPORARY VIEW IF NOT EXISTS pmid AS
SELECT DISTINCT
    genes.gene_id AS entrez,
    omim.omim_id AS omim,
    pubmed.pubmed_id AS pmid
FROM genes
JOIN omim ON genes._id = omim._id
LEFT OUTER JOIN pubmed ON genes._id = pubmed._id
