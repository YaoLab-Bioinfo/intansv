\name{readSvseq}
\alias{readSvseq}

\title{Read in the structural variations predicted by SVseq2}

\description{
    Reading the structural variations predicted by SVseq2, 
    filtering low quality predictions and merging overlapping predictions.
}

\usage{
    readSvseq(dataDir=".", regSizeLowerCutoff=100, method="SVseq2", 
              regSizeUpperCutoff=1000000, readsSupport=3) 
}

\arguments{
    \item{dataDir}{a directory containing the predictions of SVseq2.}
    \item{regSizeLowerCutoff}{the minimum size for a structural 
                              variation to be read.}
    \item{regSizeUpperCutoff}{the maximum size for a structural 
                              variation to be read.}
    \item{readsSupport}{the minimum read pair support for a 
                        structural variation to be read.}
    \item{method}{a tag to assign to the result of this function.}
}

\details{
    The predicted SVs could be further filtered by the number of reads 
    supporting the occurence of a specific SV, and the predicted size 
    of SVs to get more reliable SVs. See our paper for more details.
    The directory that specified by the parameter "dataDir" should only 
    contain the predictions of SVseq2. The deletions output files should 
    be named using the suffix ".del". See the example dataset for more details.
}

\value{
    A list with the following components:

    \item{del}{the deletions predicted by SVseq2.}
}

\author{Wen Yao}

\examples{

    svseq <- readSvseq(system.file("extdata/svseq2",package="intansv"))
    str(svseq)

}

