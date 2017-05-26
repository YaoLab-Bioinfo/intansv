\name{svAnnotation}
\alias{svAnnotation}

\title{Annotation of structural variations}

\description{
    Annotate the effect caused by structural variations to genes and 
    elements of genes.
}

\usage{
    svAnnotation(structuralVariation,genomeAnnotation)
}

\arguments{
    \item{structuralVariation}{A data frame of structural variations.}
    \item{genomeAnnotation}{A data frame of genome annotations.}
}

\details{
    A structural variation (deletion, duplication, inversion et al.) could
    affect the structure of a specific gene, including deletion of 
    introns/exons, deletion of whole gene, et al.. This function gives the 
    detailed effects caused by structural variations to genes and elements 
    of genes.

    The parameter "structuralVariation" should be a data frame with three 
    columns:
    \itemize{
         \item{chromosome}{ the chromosome of a structural variation.}
         \item{pos1}{ the start coordinate of a structural variation.}
         \item{pos2}{ the end coordinate of a structural variation.}
    }

}

\value{
    A data frame with the following columns:
    
    \item{chromosome}{the chromosome of a structural variation.}
    \item{pos1}{the start coordinate of a structural variation.}
    \item{pos2}{the end coordinate of a structural variation.}
	\item{size}{the size of a structural variation.}
	\item{info}{information on a structural variation.}
    \item{tag}{the tag of a genomic element.}
	\item{start}{the start coordinate of a genomic element.}
	\item{end}{the end coordinate of a genomic element.}
	\item{strand}{the strand of a genomic element.}
	\item{ID}{the ID of a genomic element.}
}

\author{Wen Yao}

\examples{

    breakdancer <- readBreakDancer(system.file("extdata/ZS97.breakdancer.sv",
                                   package="intansv"))
    str(breakdancer)

    msu_gff_v7 <- read.table(system.file("extdata/chr05_chr10.anno.txt", package="intansv"),
                               head=TRUE, as.is=TRUE, sep="\t")
    breakdancer.anno <- llply(breakdancer,svAnnotation,
                              genomeAnnotation=msu_gff_v7)
    str(breakdancer.anno)

}

