\name{TPGLDA}
\alias{TPGLDA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
TPGLDA(known_lncRNA_disease_interaction, known_gene_disease_interaction, gama = 0.6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{known_lncRNA_disease_interaction}{
%%     ~~Describe \code{known_lncRNA_disease_interaction} here~~
}
  \item{known_gene_disease_interaction}{
%%     ~~Describe \code{known_gene_disease_interaction} here~~
}
  \item{gama}{
%%     ~~Describe \code{gama} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (known_lncRNA_disease_interaction, known_gene_disease_interaction, 
    gama = 0.6) 
{
    A = known_lncRNA_disease_interaction
    A <- as.matrix(A)
    nl = dim(A)[1]
    nd = dim(A)[2]
    W = matrix(0, nrow = nl, ncol = nl)
    Rsocre1_lnc_dis = matrix(0, nrow = nl, ncol = nd)
    for (j in 1:nl) {
        q1 = repmat(A[j, ], nl, 1) * A/repmat(colSums(A), nl, 
            1)
        W[j, ] = 1/rowSums(A) * rowSums(q1)
    }
    Wji = t(W)
    Wj_i = Wji/repmat(colSums(Wji), nl, 1)
    Rsocre1_lnc_dis = (W + Wj_i) \%*\% A
    B = known_gene_disease_interaction
    B <- as.matrix(B)
    ng = dim(B)[1]
    W_d = matrix(0, nrow = nl, ncol = ng)
    Rsocre2_lnc_dis = matrix(0, nrow = nl, ncol = ng)
    for (j in 1:nl) {
        q2 = repmat(A[j, ], ng, 1) * B/repmat(colSums(B), ng, 
            1)
        W_d[j, ] = 1/rowSums(B) * rowSums(q2)
    }
    Rsocre2_lnc_dis = W_d \%*\% B
    final_Rscore = gama * Rsocre1_lnc_dis + (1 - gama) * Rsocre2_lnc_dis
    write.table(final_Rscore, file = "final_Rscore.txt", row.names = F, 
        quote = F, sep = "\t")
    return(final_Rscore)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
