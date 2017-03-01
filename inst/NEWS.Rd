\name{Organism.dplyr-NEWS}
\title{Organism.dplyr News}

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
