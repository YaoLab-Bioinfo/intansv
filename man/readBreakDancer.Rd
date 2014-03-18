\name{readBreakDancer}
\alias{readBreakDancer}

\title{Read in the structural variations predicted by breakDancer}

\description{
    Reading in the structural variations predicted by breakDancer, 
    filtering low quality predictions and merging overlapping predictions.
}

\usage{
    readBreakDancer(file="", scoreCutoff=60, readsSupport=3, 
                    regSizeLowerCutoff=100, regSizeUpperCutoff=1000000,
                    method="BreakDancer", ...)
}

\arguments{
    \item{file}{the output file of breakDancer.}
    \item{scoreCutoff}{the minimum score for a structural variation 
                       to be read in.}
    \item{readsSupport}{the minimum read pair support for a 
                        structural variation to be read in.}
    \item{regSizeLowerCutoff}{the minimum size for a structural 
                              variation to be read in.}
    \item{regSizeUpperCutoff}{the maximum size for a structural 
                              variation to be read in.}
    \item{method}{a tag to assign to the result of this function.}
    \item{...}{parameters passed to read.table.}
}

\details{
    The predicted SVs could be further filtered by score, 
    number of read pairs supporting the occurence of a specific SV, 
    and the predicted size of SVs to get more reliable SVs. See our 
    paper for more details.
}

\value{
    A list with the following components:

    \item{del}{the deletions predicted by breakDancer.}
    \item{inv}{the inversions predicted by breakDancer.}

}

\author{Wen Yao}

\examples{

    breakdancer <- readBreakDancer(system.file("extdata/ZS97.breakdancer.sv",
                                   package="intansv"))
    str(breakdancer)

}

