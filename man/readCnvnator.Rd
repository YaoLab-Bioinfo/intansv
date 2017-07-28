\name{readCnvnator}
\alias{readCnvnator}

\title{Read in the structural variations predicted by CNVnator}

\description{
    Reading the structural variations predicted by CNVnator, 
    filtering low quality predictions and merging overlapping predictions.
}

\usage{
    readCnvnator(dataDir=".", regSizeLowerCutoff=100, regSizeUpperCutoff=1000000,
                 method="CNVnator")
}

\arguments{
    \item{dataDir}{the directory that contain the output files of CNVnator.}
    \item{regSizeLowerCutoff}{the minimum size for a structural 
                              variation to be read.}
    \item{regSizeUpperCutoff}{the maximum size for a structural 
                              variation to be read.}
    \item{method}{a tag to assign to the result of this function.}
}

\details{
    The predicted SVs could be further filtered by the predicted 
    size of SVs to get more reliable SVs. See our paper for more details.
    The directory that specified by the parameter "dataDir" should only 
    contain the predictions of CNVnator. See the example dataset 
    for more details.
}

\value{
    A list with the following components:

    \item{del}{the deletions predicted by CNVnator.}
    \item{dup}{the duplications predicted by CNVnator.}
}

\author{Wen Yao}

\examples{

    cnvnator <- readCnvnator(system.file("extdata/cnvnator",package="intansv"))
    str(cnvnator)

}

