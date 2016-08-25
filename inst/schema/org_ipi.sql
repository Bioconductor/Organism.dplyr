CREATE TEMPORARY VIEW IF NOT EXISTS ipi AS
SELECT DISTINCT
    genes.gene_id AS entrez,
    pfam.ipi_id AS ipi,
    pfam.pfam_id AS pfam,
    prosite.prosite_id AS prosite
FROM genes
LEFT OUTER JOIN pfam ON genes._id = pfam._id
LEFT OUTER JOIN prosite ON genes._id = prosite._id
    AND pfam.ipi_id = prosite.ipi_id
