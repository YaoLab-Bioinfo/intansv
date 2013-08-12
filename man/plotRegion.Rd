\name{plotRegion}
\alias{plotRegion}

\title{Display structural variations in a specific genomic region}

\description{
    Display the structural variations in a specific genomic region 
    in circular view.
}

\usage{
    plotRegion(structuralVariation, genomeAnnotation, 
               regionChromosome, regionStart, regionEnd)
}

\arguments{
    \item{structuralVariation}{A list of structural variations.}
    \item{genomeAnnotation}{A genomic ranges of the genome annotation.}
    \item{regionChromosome}{The chromosome identifier of a specific 
                            region to view.}
    \item{regionStart}{The start coordinate of a specific region to view.}
    \item{regionEnd}{The end coordinate of a specific region to view.}
}

\details{
    Different SVs were shown as rectangles in different layers. 
    See the package vignette and the example dataset for more details.
}

\value{
    A circular plot of all the structural variations and genes in a 
    specific region with four layers:

    \itemize{
        \item{}{The composition of genes of a specific genomic region.}
        \item{}{The composition of deletions of a specific genomic region.}
        \item{}{The composition of duplications of a specific genomic region.}
        \item{}{The composition of inversions of a specific genomic region.}
    }
}

\author{Wen Yao}

\examples{

    delly <- readDelly(system.file("extdata/delly",package="intansv"))
    str(delly)

    load(system.file("extdata/genome.anno.RData",package="intansv"))
    str(msu_gff_v7)

    plotRegion(delly,msu_gff_v7,"chr05",1,200000)

}

