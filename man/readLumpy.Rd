\name{readLumpy}
\alias{readLumpy}

\title{Read in the structural variations predicted by Lumpy}

\description{
    Reading the structural variations predicted by Lumpy, 
    filtering low quality predictions and merging overlapping predictions.
}

\usage{
    readLumpy(file="", regSizeLowerCutoff=100, readsSupport=3, 
              method="Lumpy", regSizeUpperCutoff=1000000, 
              breakpointThres=200, scoreCut=0.1, ...)
}

\arguments{
    \item{file}{the file containing the prediction results of Lumpy.}
    \item{regSizeLowerCutoff}{the minimum size for a structural variation 
                              to be read.}
    \item{regSizeUpperCutoff}{the maximum size for a structural variation 
                              to be read.}
    \item{readsSupport}{the minimum read pair support for a structural variation 
                              to be read.} 
    \item{method}{a tag to assign to the result of this function.}
    \item{breakpointThres}{a threshold to remove SVs with breakpoint with 
                           too large interval.}
    \item{scoreCut}{predictions with score larger than this threshold will
                    be discarded.}
    \item{...}{parameters passed to read.table.}
}

\details{
    The predicted SVs could be further filtered by the number of reads 
    supporting the occurence of a specific SV, and the predicted size of 
    SVs to get more reliable SVs. See our paper for more details.
    The directory that specified by the parameter "dataDir" should only 
    contain the predictions of Lumpy. The deletions output files should 
    be named using the suffix "_D", the duplications output files should 
    be named using the suffix "_TD", and the inversions output files should 
    be named using the suffix "_INV". See the example dataset for more details.
}

\value{
    A list with the following components:

    \item{del}{the deletions predicted by Lumpy.}
    \item{dup}{the duplications predicted by Lumpy.}
    \item{inv}{the inversions predicted by Lumpy.}
}

\author{Wen Yao}

\examples{

    lumpy <- readLumpy(system.file("extdata/ZS97.lumpy.pesr.bedpe",package="intansv"))
    str(lumpy)

}

