
\name{geneAnnotation}
\alias{geneAnnotation}

\title{Annotation of genes affected by structural variations}

\description{
    Report the details of genes affected by structural variations.
}

\usage{
    geneAnnotation(structuralVariation,genomeAnnotation)
}

\arguments{
    \item{structuralVariation}{A data frame of structural variations.}
    \item{genomeAnnotation}{A genomic ranges of the genome annotation.}
}

\details{
    A structural variation (deletion, duplication, inversion et al.) could
    affect the structure of a specific gene, including deletion of 
    introns/exons, deletion of whole gene, et al.. And a specific gene might 
    be affected by multiple SVs. This function gives the detailed effects 
    caused by structural variations to genes and its elements from the 
    point of genes.

    The parameter "structuralVariation" should be a data frame with 
    three columns:

    \itemize{
        \item{chr}{ the chromosome of a structural variation.}
        \item{start}{ the start coordinate of a structural variation.}
        \item{end}{ the end coordinate of a structural variation.}
    }

}

\value{
    A data frame with the following columns:
    
    \item{locus}{the gene affected by structural variations.}
    \item{exon}{the effect of structural variations to exons of a specific gene.}
    \item{intron}{the effect of structural variations to introns of a 
                  specific gene.}
    \item{cds}{the effect of structural variations to cdss of a specific gene.}
    \item{utr}{the effect of structural variations to utrs of a specific gene.}

}

\author{Wen Yao}

\examples{
    breakdancer <- readBreakDancer(system.file("extdata/ZS97.breakdancer.sv",
                                   package="intansv"))
    str(breakdancer)

    load(system.file("extdata/genome.anno.RData",package="intansv"))
    str(msu_gff_v7)
    gene.breakdancer.anno <- llply(breakdancer,geneAnnotation,
                                   genomeAnnotation=msu_gff_v7)
    str(gene.breakdancer.anno)

}






