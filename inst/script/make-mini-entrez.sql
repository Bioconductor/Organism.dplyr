DELETE FROM id WHERE entrez IS NULL OR entrez NOT IN (SELECT DISTINCT entrez FROM id LIMIT 50);
SELECT COUNT(DISTINCT entrez) FROM id;
DELETE FROM id_accession WHERE entrez IS NULL OR entrez NOT IN (SELECT entrez FROM id);
DELETE FROM id_go WHERE entrez IS NULL OR entrez NOT IN (SELECT entrez FROM id);
DELETE FROM id_go_all WHERE entrez IS NULL OR entrez NOT IN (SELECT entrez FROM id);
DELETE FROM id_omim_pm WHERE entrez IS NULL OR entrez NOT IN (SELECT entrez FROM id);
DELETE FROM id_protein WHERE entrez IS NULL OR entrez NOT IN (SELECT entrez FROM id);
DELETE FROM id_transcript WHERE entrez IS NULL OR entrez NOT IN (SELECT entrez FROM id);
DELETE FROM ranges_cds WHERE entrez IS NULL OR entrez NOT IN (SELECT entrez FROM id);
DELETE FROM ranges_exon WHERE entrez IS NULL OR entrez NOT IN (SELECT entrez FROM id);
DELETE FROM ranges_gene WHERE entrez IS NULL OR entrez NOT IN (SELECT entrez FROM id);
DELETE FROM ranges_tx WHERE entrez IS NULL OR entrez NOT IN (SELECT entrez FROM id);

-- seqinfo
-- metadata_org
-- metadata_txdb

DROP index IF EXISTS cdsid_ranges_cds;
DROP index IF EXISTS ensembl_id;
DROP index IF EXISTS entrez_go;
DROP index IF EXISTS entrez_go_all;
DROP index IF EXISTS entrez_id;
DROP index IF EXISTS entrez_omim_pm;
DROP index IF EXISTS entrez_ranges_cds;
DROP index IF EXISTS entrez_ranges_exon;
DROP index IF EXISTS entrez_ranges_gene;
DROP index IF EXISTS entrez_ranges_tx;
DROP index IF EXISTS entrez_transcript;
DROP index IF EXISTS exonid_ranges_exon;
DROP index IF EXISTS symbol_id;
DROP index IF EXISTS txid_ranges_cds;
DROP index IF EXISTS txid_ranges_exon;
DROP index IF EXISTS txid_ranges_tx;

VACUUM;
