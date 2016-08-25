CREATE TEMPORARY VIEW IF NOT EXISTS view_go AS
SELECT DISTINCT
    genes.gene_id AS entrez,
    go.go_id AS go,
    go.evidence AS evidence,
    go.ontology AS ontology
FROM genes
LEFT OUTER JOIN go ON genes._id = go._id
