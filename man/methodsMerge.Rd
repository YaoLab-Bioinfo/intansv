\name{methodsMerge}
\alias{methodsMerge}

\title{Integrate structural variations predicted by different methods}

\description{
    Integrate predictions of different tools to provide more reliable 
    structural variations.
}

\usage{
    methodsMerge(..., others=NULL, 
    overLapPerDel = 0.8, overLapPerDup = 0.8, overLapPerInv = 0.8,
    numMethodsSupDel = 2, numMethodsSupDup = 2, numMethodsSupInv = 2)
}

\arguments{
    \item{...}{results of different SVs predictions read in to R by intansv.}
    \item{others}{a data frame of structural variations predicted by other tools.}
    \item{overLapPerDel}{Deletions predicted by different methods that have reciprocal coordinate overlap larger than this threshold would be clustered together}
    \item{overLapPerDup}{Duplications predicted by different methods that have reciprocal coordinate overlap larger than this threshold would be clustered together}
    \item{overLapPerInv}{Inversions predicted by different methods that have reciprocal coordinate overlap larger than this threshold would be clustered together}
    \item{numMethodsSupDel}{Deletion clusters supportted by no more than this 
    threshold of read support would be discarded}
    \item{numMethodsSupDup}{Duplication clusters supportted by no more than this
    threshold of read support would be discarded}
    \item{numMethodsSupInv}{Inversion clusters supportted by no more than this 
    threshold of read support would be discarded}
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
                                 overLapPerDel=0.7)
    str(sv_all_methods.1)

    sv_all_methods.2 <- methodsMerge(breakdancer,pindel,cnvnator,delly,svseq,
                                 overLapPerDel=0.8, numMethodsSupDel=3)
    str(sv_all_methods.2)
}

