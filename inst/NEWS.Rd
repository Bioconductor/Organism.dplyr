\name{Organism.dplyr-NEWS}
\title{Organism.dplyr News}

\section{CHANGES IN VERSION 1.32.0}{
  \subsection{BUG FIXES}{
    \itemize{
      \item (v. 1.31.1) avoid non-generic arguments in tbl.src_organism,
      closing
      \url{https://github.com/Bioconductor/Organism.dplyr/issues/19}
    }
  }
}

\section{CHANGES IN VERSION 1.30.0}{
  \subsection{BUG FIXES}{
    \itemize{
      \item (v. 1.30.1) avoid non-generic arguments in tbl.src_organism,
      closing
      \url{https://github.com/Bioconductor/Organism.dplyr/issues/19}

      \item \code{src_ucsc()} failed to correctly handle new 'hs1'
      resources (T2T genomes) for 'human'
    }
  }
}

\section{CHANGES IN VERSION 1.14.0}{
  \subsection{NEW FEATURES}{
    \itemize{

      \item \code{src_organism()} supports an option
      \code{overwrite=FALSE} to optionally over-write exisiting (cached)
      resources created from a previous txdb version.

      \item \code{src_organism()} supports construction from a TxDb object.

    }
  }
}

\section{CHANGES IN VERSION 1.0.0}{

  \subsection{NEW FEATURES}{
    \itemize{

      \item This package provides an interface to combined
      \emph{Bioconductor} \code{org.*} (identifier) and \code{TxDb.*}
      (genomic coordinate) annotation resources. The interface is
      implemented at several levels, including low-level \sQuote{dplyr},
      org-like \code{select()}, and TxDb-like \code{genes()}, etc.

    }
  }

  \subsection{SIGNIFICANT USER-VISIBLE CHANGES}{
    \itemize{

      \item Filters use strict CamelCase convention.

    }
  }

  \subsection{BUG FIXES}{
    \itemize{

      \item \code{*IdFilter} and \code{*RankFilter} are numeric
      (integer), rather than character.

    }
  }

}
