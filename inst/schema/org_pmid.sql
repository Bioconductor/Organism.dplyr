CREATE TEMPORARY VIEW IF NOT EXISTS pmid AS
SELECT DISTINCT
    genes.gene_id AS entrez,
    pubmed.pubmed_id AS pmid
FROM genes
LEFT OUTER JOIN pubmed ON genes._id = pubmed._id
