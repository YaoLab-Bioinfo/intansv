\name{readDelly}
\alias{readDelly}

\title{Read in the structural variations predicted by DELLY}

\description{
    Reading the structural variations predicted by DELLY, 
    filtering low quality predictions and merging overlapping predictions.
}

\usage{
    readDelly(dataDir=".", regSizeLowerCutoff=100, 
              regSizeUpperCutoff=1000000, readsSupport=3,
              method="DELLY", pass=TRUE, minMappingQuality=30)
}

\arguments{
    \item{dataDir}{a directory containing the prediction results of DELLY.}
    \item{regSizeLowerCutoff}{the minimum size for a structural 
                              variation to be read.}
    \item{regSizeUpperCutoff}{the maximum size for a structural 
                              variation to be read.}
    \item{readsSupport}{the minimum read pair support for a structural 
                        variation to be read.}
    \item{method}{a tag to assign to the result of this function.}
    \item{pass}{set pass=TRUE to remove LowQual SVs reported by DELLY.}
    \item{minMappingQuality}{the minimum mapping quality for a SV to be read.}
}

\details{
    The predicted SVs could be further filtered by the number of read pairs 
    supporting the occurence of a specific SV, and the predicted size of SVs 
    to get more reliable SVs. See our paper for more details.
    The directory that specified by the parameter "dataDir" should only 
    contain the predictions of DELLY. See the example dataset for more details.
}

\value{
    A list with the following components:

    \item{del}{the deletions predicted by DELLY.}
    \item{dup}{the duplications predicted by DELLY.}
    \item{inv}{the inversions predicted by DELLY.}
}

\author{Wen Yao}

\examples{

    delly <- readDelly(system.file("extdata/delly",package="intansv"))
    str(delly)

}

