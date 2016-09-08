CREATE TEMPORARY VIEW IF NOT EXISTS id_accession AS
SELECT DISTINCT
    genes.gene_id AS entrez,
    accessions.accession AS accnum,
    refseq.accession AS refseq
FROM genes
LEFT OUTER JOIN accessions ON genes._id = accessions._id
LEFT OUTER JOIN refseq ON genes._id = refseq._id
    AND refseq.accession = accessions.accession;

CREATE TEMPORARY VIEW IF NOT EXISTS id_transcript AS
SELECT DISTINCT
    genes.gene_id AS entrez,
    unigene.unigene_id AS unigene,
    ensembl_trans.trans_id AS ensembltrans
FROM genes
LEFT OUTER JOIN unigene ON genes._id = unigene._id
LEFT OUTER JOIN ensembl_trans ON genes._id = ensembl_trans._id;

CREATE TEMPORARY VIEW IF NOT EXISTS id AS
SELECT DISTINCT
    genes.gene_id AS entrez,
    cytogenetic_location AS map,
    ensembl.ensembl_id AS ensembl,
    gene_info.symbol AS symbol,
    gene_info.gene_name AS genename,
    alias.alias_symbol AS alias
FROM genes
JOIN gene_info ON genes._id = gene_info._id
JOIN cytogenetic_locations ON genes._id = cytogenetic_locations._id
JOIN ensembl ON genes._id = ensembl._id
LEFT OUTER JOIN alias ON genes._id = alias._id;

CREATE TEMPORARY VIEW IF NOT EXISTS id_omim_pm AS
SELECT DISTINCT
    genes.gene_id AS entrez,
    omim.omim_id AS omim,
    pubmed.pubmed_id AS pmid
FROM genes
JOIN omim ON genes._id = omim._id
LEFT OUTER JOIN pubmed ON genes._id = pubmed._id;

CREATE TEMPORARY VIEW IF NOT EXISTS id_protein AS
SELECT DISTINCT
    genes.gene_id AS entrez,
    ec.ec_number AS enzyme,
    ensembl_prot.prot_id AS ensemblprot,
    uniprot.uniprot_id AS uniprot,
    pfam.ipi_id AS ipi,
    pfam.pfam_id AS pfam,
    prosite.prosite_id AS prosite
FROM genes
LEFT OUTER JOIN ec ON genes._id = ec._id
LEFT OUTER JOIN uniprot ON genes._id = uniprot._id
LEFT OUTER JOIN ensembl_prot ON genes._id = ensembl_prot._id
LEFT OUTER JOIN pfam ON genes._id = pfam._id
LEFT OUTER JOIN prosite ON genes._id = prosite._id
    AND pfam.ipi_id = prosite.ipi_id;

CREATE TEMPORARY VIEW IF NOT EXISTS id_go AS
SELECT DISTINCT
    main.genes.gene_id AS entrez,
    main.go.go_id AS go,
    main.go.evidence AS evidence,
    main.go.ontology AS ontology
FROM main.genes
LEFT OUTER JOIN go ON main.genes._id = main.go._id;

CREATE TEMPORARY VIEW IF NOT EXISTS id_go_all AS
SELECT DISTINCT
    main.genes.gene_id AS entrez,
    main.go_all.go_id AS goall,
    main.go_all.evidence AS evidenceall,
    main.go_all.ontology AS ontologyall
FROM main.genes
LEFT OUTER JOIN main.go_all ON main.genes._id = main.go_all._id;
