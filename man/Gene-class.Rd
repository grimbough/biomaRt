\name{Gene-class}
\doctype{class}
\alias{Gene-class}
\title{Class Gene}
\description{Represents a Gene class}
\section{Creating Objects}{
\code{ new('Gene',}\cr
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
	This Class is used to represent a Gene.  It stores the query identifier, the mart specific identifier,  gene symbol,
	a description, the chromosome and band on which the gene is located,
	the start and end positions of the gene on the chromosome.  The Gene class 
	also has slots for GO and OMIM annotation.
}
\author{Steffen Durinck, \url{http://www.esat.kuleuven.ac.be/~sdurinck}} 

\seealso{}

\examples{


}
\keyword{methods}

