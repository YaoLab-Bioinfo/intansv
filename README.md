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
<li><a href="https://github.com/abyzovlab/CNVnator" target="_blank">CNVnator</a></li>
<li><a href="https://github.com/arq5x/lumpy-sv" target="_blank">Lumpy</a></li>
<li><a href="https://github.com/Steven-N-Hart/softsearch" target="_blank">SoftSearch</a></li>
</ul>

### Install

-   the latest development version from Github with

    ``` r
    install.packages("devtools")  
    devtools::install_github("venyao/intansv")    
    ```
	
-   the latest development version from Bioconductor with

    ``` r
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

    # The following initializes usage of Bioc devel
    BiocManager::install(version='devel')
    BiocManager::install("intansv")  
    ```
	
Help Manual
------------
[Overview of intansv](https://github.com/venyao/intansv/blob/master/intansvOverview.pdf)
    
Example output of these SV programs are listed at https://venyao.github.io/intansv/  
	
### intansv development version 1.27.0 support the following version of SV prediction programs
* BreakDancer version 1.4.5  
* Pindel version v0.2.5b9   
* DELLY version 0.8.1  
* SVseq2 version 2.2  
* CNVnator version 0.3.3  
* Lumpy version 0.3.0  
* SoftSearch version 2.4  

### intansv development version 1.21.0 support the following version of SV prediction programs
* BreakDancer version 1.4.5  
* Pindel version v0.2.5b9   
* DELLY version 0.7.6  
* SVseq2 version 2.2  
* CNVnator version 0.3.3  
* Lumpy version 0.2.13  
* SoftSearch version 2.4  

### intansv release version 1.7.3 support the following version of SV prediction programs  
* BreakDancer version 1.4.5
* Pindel version 0.2.5a8
* DELLY version 0.6.1
* SVseq2 version 2.2
* CNVnator version 0.3
* Lumpy version 0.2.8
* SoftSearch version 2.4

### intansv release version 1.7.1 support the following version of SV prediction programs
* BreakDancer version 1.0r112
* Pindel version 0.2.4t
* DELLY version 0.0.11
* SVseq2 version 2.2
* CNVnator version 0.2.7
* Lumpy version 0.1.5
* SoftSearch version 2.3

Please go to [Bioconductor](http://www.bioconductor.org) to download the old versions of the intansv package.  
