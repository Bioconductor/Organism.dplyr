% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter.R
\docType{class}
\name{BasicFilter-class}
\alias{BasicFilter-class}
\alias{AccnumFilter}
\alias{AliasFilter}
\alias{CdsChromFilter}
\alias{CdsIdFilter}
\alias{CdsNameFilter}
\alias{CdsStrandFilter}
\alias{EnsemblFilter}
\alias{EnsemblprotFilter}
\alias{EnsembltransFilter}
\alias{EnzymeFilter}
\alias{EvidenceFilter}
\alias{EvidenceallFilter}
\alias{ExonChromFilter}
\alias{ExonStrandFilter}
\alias{FlybaseFilter}
\alias{FlybaseCgFilter}
\alias{FlybaseProtFilter}
\alias{GeneChromFilter}
\alias{GeneStrandFilter}
\alias{GoFilter}
\alias{GoallFilter}
\alias{IpiFilter}
\alias{MapFilter}
\alias{MgiFilter}
\alias{OmimFilter}
\alias{OntologyFilter}
\alias{OntologyallFilter}
\alias{PfamFilter}
\alias{PmidFilter}
\alias{PrositeFilter}
\alias{RefseqFilter}
\alias{TxChromFilter}
\alias{TxStrandFilter}
\alias{TxTypeFilter}
\alias{WormbaseFilter}
\alias{ZfinFilter}
\alias{CharacterFilter-class}
\alias{IntegerFilter-class}
\alias{show,CharacterFilter-method}
\alias{show,IntegerFilter-method}
\alias{show,BasicFilter-method}
\alias{supportedFilters,src_organism-method}
\title{Filtering src_organism objects}
\usage{
AccnumFilter(value, condition = "==")
AliasFilter(value, condition = "==")
CdsChromFilter(value, condition = "==")
CdsIdFilter(value, condition = "==")
CdsNameFilter(value, condition = "==")
CdsStrandFilter(value, condition = "==")
EnsemblFilter(value, condition = "==")
EnsemblprotFilter(value, condition = "==")
EnsembltransFilter(value, condition = "==")
EnzymeFilter(value, condition = "==")
EvidenceFilter(value, condition = "==")
EvidenceallFilter(value, condition = "==")
ExonChromFilter(value, condition = "==")
ExonStrandFilter(value, condition = "==")
FlybaseFilter(value, condition = "==")
FlybaseCgFilter(value, condition = "==")
FlybaseProtFilter(value, condition = "==")
GeneChromFilter(value, condition = "==")
GeneStrandFilter(value, condition = "==")
GoFilter(value, condition = "==")
GoallFilter(value, condition = "==")
IpiFilter(value, condition = "==")
MapFilter(value, condition = "==")
MgiFilter(value, condition = "==")
OmimFilter(value, condition = "==")
OntologyFilter(value, condition = "==")
OntologyallFilter(value, condition = "==")
PfamFilter(value, condition = "==")
PmidFilter(value, condition = "==")
PrositeFilter(value, condition = "==")
RefseqFilter(value, condition = "==")
TxChromFilter(value, condition = "==")
TxStrandFilter(value, condition = "==")
TxTypeFilter(value, condition = "==")
WormbaseFilter(value, condition = "==")
ZfinFilter(value, condition = "==")

\S4method{show}{BasicFilter}(object)

\S4method{supportedFilters}{src_organism}(object)
}
\arguments{
\item{object}{A \code{BasicFilter} or \code{GRangesFilter} object}

\item{value}{Value of the filter. For \code{GRangesFilter} value should be a
\code{GRanges} object.}

\item{condition}{The condition to be used in filter for genomic
extractors, one of "==", "!=", "startsWith", "endsWith", ">",
"<", ">=", "<=".  For character values "==", "!=", "startsWith"
and "endsWith" are allowed, for numeric values
(\code{CdsStartFilter}, \code{CdsEndFilter},
\code{ExonStartFilter}, \code{ExonEndFilter},
\code{GeneStartFilter}, \code{GeneEndFilter},
\code{TxStartFilter} and \code{TxEndFilter}), "==", "!=", ">",
">=", "<" and "<=". Default condition is "==".}
}
\value{
A Filter object showing class, value and condition of the filter
}
\description{
These functions create filters to be used by the "select"
interface to \code{src_organism} objects.
}
\details{
All filters except \code{GRangesFilter()} takes value(s) from
corresponding fields in the data base. For example,
\code{AccnumFilter()} takes values of accession number(s), which
come from field \code{accnum}. See \code{keytypes()} and
\code{keys()} for possible values.

\code{GRangesFilter()} takes a \code{GRanges} object as filter, and returns
genomic extractors (\code{genes}, \code{transcripts}, etc.) that are
partially overlapping with the region.

\code{supportedFilters()} lists all available filters for
\code{src_organism} object.
}
\examples{
src <- src_organism(dbpath=hg38light())
keytypes(src)
head(keys(src, "ensembl"))

## filter by ensembl
EnsemblFilter("ENSG00000171862")

## filter by gene symbol start with "BRAC"
SymbolFilter("BRCA", "startsWith")

## filter by GRanges
GRangesFilter(GenomicRanges::GRanges("chr10:87869000-87876000"))

## filter by transcript start position
TxStartFilter(87863438, ">")

}
\seealso{
\code{\link{src_organism}} for creating a \code{src_organism}
    object.

    \code{\link{transcripts_tbl}} for generic functions to extract
     genomic features from a \code{src_organism} object.

     \code{\link{select,src_organism-method}} for "select"
    interface on \code{src_organism} objects.
}
\author{
Yubo Cheng.
}
