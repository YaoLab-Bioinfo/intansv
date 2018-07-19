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
<<<<<<< HEAD
<li><a href="https://github.com/abyzovlab/CNVnator" target="_blank">CNVnator</a></li>
<li><a href="https://github.com/arq5x/lumpy-sv" target="_blank">Lumpy</a></li>
<li><a href="https://github.com/Steven-N-Hart/softsearch" target="_blank">SoftSearch</a></li>
=======
<li><a href="http://sv.gersteinlab.org/cnvnator/" target="_blank">CNVnator</a></li>
<li><a href="https://github.com/arq5x/lumpy-sv" target="_blank">Lumpy</a></li>
<li><a href="http://code.google.com/p/softsearch/" target="_blank">SoftSearch</a></li>
>>>>>>> upstream/master
</ul>

### Install

<<<<<<< HEAD
-   the latest development version from Github with

    ``` r
    install.packages("devtools")  
    devtools::install_github("venyao/intansv")    
    ```
	
-   the latest development version from Bioconductor with
=======
-   the latest released version from Bioconductor with

    ``` r
    source("http://bioconductor.org/biocLite.R")
    biocLite("intansv")
    ```
    
-   the latest development version from github with
>>>>>>> upstream/master

    ``` r
    source("http://bioconductor.org/biocLite.R")  
    useDevel()  
    biocLite("intansv")  
    ```
<<<<<<< HEAD
	
Help Manual
------------
[Overview of intansv](https://github.com/venyao/intansv/blob/master/intansvOverview.pdf)
    
Example output of these SV programs are listed at http://venyao.github.io/intansv/  
	
### intansv development version 1.17.1 support the following version of SV prediction programs
* BreakDancer version 1.4.5  
* Pindel version v0.2.5b9   
* DELLY version 0.7.6  
* SVseq2 version 2.2  
* CNVnator version 0.3.3  
* Lumpy version 0.2.13  
* SoftSearch version 2.4  
=======
>>>>>>> upstream/master

### intansv release version 1.7.3 support the following version of SV prediction programs  
* BreakDancer version 1.4.5
* Pindel version 0.2.5a8
* DELLY version 0.6.1
* SVseq2 version 2.2
* CNVnator version 0.3
* Lumpy version 0.2.8
* SoftSearch version 2.4

<<<<<<< HEAD
=======
Example output of these SV programs are listed at http://andrewhzau.github.io/intansv/  

>>>>>>> upstream/master
### intansv release version 1.7.1 support the following version of SV prediction programs
* BreakDancer version 1.0r112
* Pindel version 0.2.4t
* DELLY version 0.0.11
* SVseq2 version 2.2
* CNVnator version 0.2.7
* Lumpy version 0.1.5
* SoftSearch version 2.3


