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
    \item{genomeAnnotation}{A genomic ranges of the genome annotation.}
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
    
    \item{chr}{the chromosome of a structural variation.}
    \item{start}{the start coordinate of a structural variation.}
    \item{end}{the end coordinate of a structural variation.}
    \item{overlap}{the overlap length between a structural variation and a 
                   specific gene or its element.}
    \item{annotation}{the annotation of a specific gene that overlap with 
                      the structural variation.}
    \item{parent}{the ID of a specific gene that overlap with the 
                  structural variation.}

}

\author{Wen Yao}

\examples{

    breakdancer <- readBreakDancer(system.file("extdata/ZS97.breakdancer.sv",
                                   package="intansv"))
    str(breakdancer)

    load(system.file("extdata/genome.anno.RData",package="intansv"))
    str(msu_gff_v7)
    breakdancer.anno <- llply(breakdancer,svAnnotation,
                              genomeAnnotation=msu_gff_v7)
    str(breakdancer.anno)

}

