CREATE TEMPORARY VIEW IF NOT EXISTS view_go_all AS
SELECT DISTINCT
    genes.gene_id AS entrez,
    go_all.go_id AS goall,
    go_all.evidence AS evidenceall,
    go_all.ontology AS ontologyall
FROM genes
LEFT OUTER JOIN go_all ON genes._id = go_all._id
