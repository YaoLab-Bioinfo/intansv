\name{readPindel}
\alias{readPindel}

\title{Read in the structural variations predicted by Pindel}

\description{
    Reading the structural variations predicted by Pindel, 
    filtering low quality predictions and merging overlapping predictions.
}

\usage{
    readPindel(dataDir=".", regSizeLowerCutoff=100, 
               regSizeUpperCutoff=1000000, readsSupport=3,
               method="Pindel")
}

\arguments{
    \item{dataDir}{the directory containing the prediction results of Pindel.}
    \item{regSizeLowerCutoff}{the minimum size for a structural variation 
                              to be read.}
    \item{regSizeUpperCutoff}{the maximum size for a structural variation 
                              to be read.}
    \item{readsSupport}{the minimum read pair support for a structural 
                        variation to be read.}
    \item{method}{a tag to assign to the result of this function.}
}

\details{
    The predicted SVs could be further filtered by the number of reads 
    supporting the occurence of a specific SV, and the predicted size of 
    SVs to get more reliable SVs. See our paper for more details.
    The directory that specified by the parameter "dataDir" should only 
    contain the predictions of Pindel. The deletions output files should 
    be named using the suffix "_D", the duplications output files should 
    be named using the suffix "_TD", and the inversions output files should 
    be named using the suffix "_INV". See the example dataset for more details.
}

\value{
    A list with the following components:

    \item{del}{the deletions predicted by Pindel.}
    \item{dup}{the duplications predicted by Pindel.}
    \item{inv}{the inversions predicted by Pindel.}
}

\author{Wen Yao}

\examples{

    pindel <- readPindel(system.file("extdata/pindel",package="intansv"))
    str(pindel)

}

