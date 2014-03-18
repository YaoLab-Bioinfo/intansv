\name{readSoftSearch}
\alias{readSoftSearch}

\title{Read in the structural variations predicted by SoftSearch}

\description{
    Reading the structural variations predicted by SoftSearch, 
    filtering low quality predictions and merging overlapping predictions.
}

\usage{
    readSoftSearch(file="", regSizeLowerCutoff=100, readsSupport=3, 
                   method="softSearch", regSizeUpperCutoff=1000000, 
                   softClipsSupport=3, ...) 
}

\arguments{
    \item{file}{the file containing the prediction results of SoftSearch.}
    \item{regSizeLowerCutoff}{the minimum size for a structural variation 
                              to be read.}
    \item{regSizeUpperCutoff}{the maximum size for a structural variation 
                              to be read.}
    \item{readsSupport}{the minimum read pair support for a structural 
                        variation to be read.}
    \item{method}{a tag to assign to the result of this function.}
    \item{softClipsSupport}{the minimum soft clip support for a structural 
                            variation to be read.}
    \item{...}{parameters passed to read.table}
}

\details{
    The predicted SVs could be further filtered by the number of reads 
    supporting the occurence of a specific SV, and the predicted size of 
    SVs to get more reliable SVs. See our paper for more details.
    The directory that specified by the parameter "dataDir" should only 
    contain the predictions of SoftSearch. The deletions output files should 
    be named using the suffix "_D", the duplications output files should 
    be named using the suffix "_TD", and the inversions output files should 
    be named using the suffix "_INV". See the example dataset for more details.
}

\value{
    A list with the following components:

    \item{del}{the deletions predicted by SoftSearch.}
    \item{dup}{the duplications predicted by SoftSearch.}
    \item{inv}{the inversions predicted by SoftSearch.}
}

\author{Wen Yao}

\examples{

    softSearch <- readSoftSearch(system.file("extdata/ZS97.softsearch",package="intansv"))
    str(softSearch)

}

