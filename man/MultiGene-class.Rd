\name{MultiGene-class}
\doctype{class}
\alias{MultiGene-class}
\title{Class MultiGene}
\description{Represents a set of Genes}
\section{Creating Objects}{
\code{ new('MultiGene',}\cr
\code{	id = ..../, #query identifier of class character}\cr
\code{	martID = ..../, #mart specific identifier of class character}\cr
\code{	symbol = ..../, #gene symbol of class character}\cr
\code{	description = ..../, #gene description of class character}\cr
\code{	chromosome = ..../, #chromosome of class character}\cr
\code{	band = ..../, #band of class character}\cr
\code{	start = ..../, #start on chromosome of class numeric}\cr
\code{	end = ..../, #end on chromosome of class numeric}\cr
\code{	GO = ..../, #GO annotation of class GO}\cr
\code{  OMIM = ..../, #OMIM annotation of class OMIM}\cr
\code{ )}}
\details{
	This Class is used to represent a MultiGene.  It stores the query identifier, the mart specific identifiers, gene symbols,
	a descriptions, the chromosomes and bands on which the genes are located,
	the start and end positions of the genes on the chromosomes.  The MultiGene class 
	also has slots for GO and OMIM annotation.
}
\author{Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\seealso{}

\examples{


}
\keyword{methods}

