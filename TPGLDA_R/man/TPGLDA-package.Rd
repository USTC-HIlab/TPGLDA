\name{TPGLDA-package}
\alias{TPGLDA-package}
\alias{TPGLDA}
\docType{package}
\title{TPGLDA:Novel prediction of associations between lncRNAs and diseases via lncRNA-disease-gene tripartite graph
\packageTitle{TPGLDA}
}
\description{This package is implemented TPGLDA algoritm
\packageDescription{TPGLDA}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{TPGLDA}
\packageIndices{TPGLDA}
~~ input following code to run TPGLDA and achieve related predictive results ~~
~~ TPGLDA(known_lncRNA_disease_interaction,known_gene_disease_interaction,gama=0.6) ~~
}
\author{
\packageAuthor{TPGLDA}
Liang,Ding
Maintainer: \packageMaintainer{TPGLDA}
}

~~ Optionally other standard keywords, one per line, from file KEYWORDS ~~
~~ in the R documentation directory ~~
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
    known_gene_disease_interaction <- read.delim("~/data/known_gene_disease_interaction.txt", header=FALSE)
	known_lncRNA_disease_interaction <- read.delim("~/data/known_lncRNA_disease_interaction.txt", header=FALSE)
	lncRNAsimilarity <- read.delim("~/data/lncRNAsimilarity.txt", header=FALSE)
	diseasesimilarith <- read.delim("~/data/diseasesimilarith.txt", header=FALSE)
    TPGLDA(known_lncRNA_disease_interaction,known_gene_disease_interaction,gama=0.6);
}
