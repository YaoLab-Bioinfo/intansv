\name{methodsMerge}
\alias{methodsMerge}

\title{Integrate structural variations predicted by different methods}

\description{
    Integrate predictions of different tools to provide more reliable 
    structural variations.
}

\usage{
    methodsMerge(breakdancer, pindel, cnvnator, delly, svseq, others=NULL)
}

\arguments{
    \item{breakdancer}{a list of structural variations predicted by BreakDancer.}
    \item{pindel}{a list of structural variations predicted by Pindel.}
    \item{cnvnator}{a list of structural variations predicted by CNVnator.}
    \item{delly}{a list of structural variations predicted by DELLY.}
    \item{svseq}{a list of structural variations predicted by SVseq2.}
    \item{others}{a data frame of structural variations predicted by other tools.}
    \item{overLapPer}{SVs predicted by different methods that have reciprocal coordinate overlap larger than this threshold would be clustered together}
    \item{numMethodsSup}{SV clusters supportted by no more than this threshold would be discarded}
}

\details{
    A structural variation (deletion, duplication, inversion et al.) may
    be reported by different tools. However, the boundaries of this structural
    variation predicted by different tools don't always agree with each other.
    Predictions of different methods with reciprocal overlap more than 80 
    percent were merged. Structural varions supported by only one method were 
    discarded. 
}

\value{
   A list with the following components:

   \item{del}{the integrated deletions of different methods.}
   \item{dup}{the integrated duplications of different methods.}
   \item{inv}{the integrated inversions of different methods.}
}

\author{Wen Yao}

\examples{
    breakdancer <- readBreakDancer(system.file("extdata/ZS97.breakdancer.sv",
                                   package="intansv"))
    str(breakdancer)

    cnvnator <- readCnvnator(system.file("extdata/cnvnator",package="intansv"))
    str(cnvnator)

    svseq <- readSvseq(system.file("extdata/svseq2",package="intansv"))
    str(svseq)

    delly <- readDelly(system.file("extdata/delly",package="intansv"))
    str(delly)

    pindel <- readPindel(system.file("extdata/pindel",package="intansv"))
    str(pindel)

    sv_all_methods <- methodsMerge(breakdancer,pindel,cnvnator,delly,svseq)
    str(sv_all_methods)
    
    sv_all_methods.1 <- methodsMerge(breakdancer,pindel,cnvnator,delly,svseq,
                                 overLapPer=0.7)
    str(sv_all_methods.1)

    sv_all_methods.2 <- methodsMerge(breakdancer,pindel,cnvnator,delly,svseq,
                                 overLapPer=0.8, numMethodsSup=3)
    str(sv_all_methods.2)
}

