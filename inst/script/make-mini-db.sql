DELETE FROM id WHERE ensembl IS NULL OR ensembl NOT IN (SELECT DISTINCT ensembl FROM id LIMIT 50);
DELETE FROM id_accession WHERE ensembl IS NULL OR ensembl NOT IN (SELECT ensembl FROM id);
DELETE FROM id_go WHERE ensembl IS NULL OR ensembl NOT IN (SELECT ensembl FROM id);
DELETE FROM id_go_all WHERE ensembl IS NULL OR ensembl NOT IN (SELECT ensembl FROM id);
DELETE FROM id_pm WHERE ensembl IS NULL OR ensembl NOT IN (SELECT ensembl FROM id);
DELETE FROM id_protein WHERE ensembl IS NULL OR ensembl NOT IN (SELECT ensembl FROM id);
DELETE FROM id_transcript WHERE ensembl IS NULL OR ensembl NOT IN (SELECT ensembl FROM id);
DELETE FROM ranges_cds WHERE ensembl IS NULL OR ensembl NOT IN (SELECT ensembl FROM id);
DELETE FROM ranges_exon WHERE ensembl IS NULL OR ensembl NOT IN (SELECT ensembl FROM id);
DELETE FROM ranges_gene WHERE ensembl IS NULL OR ensembl NOT IN (SELECT ensembl FROM id);
DELETE FROM ranges_tx WHERE ensembl IS NULL OR ensembl NOT IN (SELECT ensembl FROM id);

-- seqinfo
-- metadata_org
-- metadata_txdb

DROP index IF EXISTS cdsid_ranges_cds;
DROP index IF EXISTS ensembl_accession;
DROP index IF EXISTS ensembl_go;
DROP index IF EXISTS ensembl_go_all;
DROP index IF EXISTS ensembl_id;
DROP index IF EXISTS ensembl_pm;
DROP index IF EXISTS ensembl_protein;
DROP index IF EXISTS ensembl_ranges_cds;
DROP index IF EXISTS ensembl_ranges_exon;
DROP index IF EXISTS ensembl_ranges_gene;
DROP index IF EXISTS ensembl_ranges_tx;
DROP index IF EXISTS ensembl_transcript;
DROP index IF EXISTS entrez_accession;
DROP index IF EXISTS entrez_go;
DROP index IF EXISTS entrez_go_all;
DROP index IF EXISTS entrez_id;
DROP index IF EXISTS entrez_pm;
DROP index IF EXISTS entrez_protein;
DROP index IF EXISTS entrez_transcript;
DROP index IF EXISTS exonid_ranges_exon;
DROP index IF EXISTS mgi_id;
DROP index IF EXISTS symbol_id;
DROP index IF EXISTS txid_ranges_cds;
DROP index IF EXISTS txid_ranges_exon;
DROP index IF EXISTS txid_ranges_tx;

VACUUM;
