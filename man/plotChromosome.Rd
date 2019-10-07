\name{plotChromosome}
\alias{plotChromosome}

\title{Display the chromosome distribution of structural variations}

\description{
    Display the chromosome distribution of structural variations by splitting 
    the chromosomes into windows of specific size and counting the number of 
    structural variations in each window.
}

\usage{
    plotChromosome(genome, structuralVariation, windowSize=1000000)
}

\arguments{
    \item{genome}{A data frame with ID and length of all Chromosomes.}
    \item{structuralVariation}{A list of structural variations.}
    \item{windowSize}{A specific size (in base pair) to split chromosomes 
                      into windows.}
}

\details{
    To visualize the distribution of structural variations in the whole genome, 
    chromosomes were splitted into windows of specific size (default 1 Mb) and 
    the number of structural variations in each window were counted. The number 
    of structural variations were shown using circular barplot.
}

\value{
   A circular plot with five layers:

\itemize{
   \item{}{the circular view of genome ideogram.}
   \item{}{the chromosome coordinates labels.}
   \item{}{the circular barplot of number of deletions in each chromosome 
           window.}
   \item{}{the circular barplot of number of duplications in each chromosome 
           window.}
   \item{}{the circular barplot of number of inversions in each chromosome 
           window.}
}
}

\author{Wen Yao}

\examples{

    delly <- readDelly(system.file("extdata/ZS97.DELLY.vcf",package="intansv"))
    str(delly)

    genome.file.path <- system.file("extdata/chr05_chr10.genome.txt", package="intansv")
    genome <- read.table(genome.file.path, head=TRUE, as.is=TRUE)
    str(genome)

    plotChromosome(genome,delly,1000000)
	
}

