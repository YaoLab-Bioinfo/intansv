intansv
=======

Integrative analysis of structural variations (SV)

This package provides efficient tools to read and integrate structural variations predicted by popular softwares. Annotation and visualization of structural variations are also implemented in the package.

### intansv integrates the SV predictions of the following programs
<ul>
<li><a href="https://github.com/genome/breakdancer" target="_blank">BreakDancer</a></li>
<li><a href="https://github.com/genome/pindel" target="_blank">Pindel</a></li>
<li><a href="https://github.com/tobiasrausch/delly" target="_blank">DELLY</a></li>
<li><a href="http://www.engr.uconn.edu/~jiz08001/svseq2.html" target="_blank">SVseq2</a></li>
<li><a href="http://sv.gersteinlab.org/cnvnator/" target="_blank">CNVnator</a></li>
<li><a href="https://github.com/arq5x/lumpy-sv" target="_blank">Lumpy</a></li>
<li><a href="http://code.google.com/p/softsearch/" target="_blank">SoftSearch</a></li>
</ul>

### Install

-   the latest released version from Bioconductor with

    ``` r
    source("http://bioconductor.org/biocLite.R")
    biocLite("intansv")
    ```
    
-   the latest development version from github with

    ``` r
    if (packageVersion("devtools") < 1.6) {
      install.packages("devtools")
    }
    devtools::install_github("andrewhzau/intansv")
    ```

### intansv release version 1.7.3 support the following version of SV prediction programs  
* BreakDancer version 1.4.5
* Pindel version 0.2.5a8
* DELLY version 0.6.1
* SVseq2 version 2.2
* CNVnator version 0.3
* Lumpy version 0.2.8
* SoftSearch version 2.4

Example output of these SV programs are listed here: http://andrewhzau.github.io/intansv/  

### intansv release version 1.7.1 support the following version of SV prediction programs
* BreakDancer version 1.0r112
* Pindel version 0.2.4t
* DELLY version 0.0.11
* SVseq2 version 2.2
* CNVnator version 0.2.7
* Lumpy version 0.1.5
* SoftSearch version 2.3


