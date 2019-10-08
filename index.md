---
layout: page
title: Welcome to intansv!
tagline: 
---
{% include JB/setup %}

__Integrative analysis of structural variations__

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

    ```r
    install.packages("devtools")  
    devtools::install_github("venyao/intansv")  
	```
	
-   the latest development version from Bioconductor with

	```r
    if (!requireNamespace("BiocManager", quietly = TRUE))  
        install.packages("BiocManager")  
    BiocManager::install(version='devel')  
    BiocManager::install("intansv")  
	```
    

### Help Manual
[Overview of intansv](https://github.com/venyao/intansv/blob/master/intansvOverview.pdf)  


### intansv development version 1.27.0 support the following version of SV prediction programs
* BreakDancer version 1.4.5, [example output]({{ BASE_PATH }}/assets/BreakDancer_version_1.4.5.output)  
* Pindel version v0.2.5b9, [example output]({{ BASE_PATH }}/assets/Pindel_version_0.2.5b9.output.zip)   
* DELLY version 0.8.1, [example output]({{ BASE_PATH }}/assets/DELLY_version_0.8.1.output)  
* SVseq2 version 2.2, [example output]({{ BASE_PATH }}/assets/SVseq2_version_2.2.output.zip)  
* CNVnator version 0.3.3, [example output]({{ BASE_PATH }}/assets/CNVnator_version_0.3.3.output.zip)  
* Lumpy version 0.3.0, [example output]({{ BASE_PATH }}/assets/Lumpy_version_0.3.0.output)    
* SoftSearch version 2.4, [example output]({{ BASE_PATH }}/assets/SoftSearch_version_2.4.output)  

### intansv development version 1.21.0 support the following version of SV prediction programs
* BreakDancer version 1.4.5, [example output]({{ BASE_PATH }}/assets/BreakDancer_version_1.4.5.output)  
* Pindel version v0.2.5b9, [example output]({{ BASE_PATH }}/assets/Pindel_version_0.2.5b9.output.zip)   
* DELLY version 0.7.6, [example output]({{ BASE_PATH }}/assets/DELLY_version_0.7.6.output.zip)  
* SVseq2 version 2.2, [example output]({{ BASE_PATH }}/assets/SVseq2_version_2.2.output.zip)  
* CNVnator version 0.3.3, [example output]({{ BASE_PATH }}/assets/CNVnator_version_0.3.3.output.zip)  
* Lumpy version 0.2.13, [example output]({{ BASE_PATH }}/assets/Lumpy_version_0.2.13.output)  
* SoftSearch version 2.4, [example output]({{ BASE_PATH }}/assets/SoftSearch_version_2.4.output)  

### intansv development version 1.17.1 support the following version of SV prediction programs
* BreakDancer version 1.4.5, [example output]({{ BASE_PATH }}/assets/BreakDancer_version_1.4.5.output)  
* Pindel version v0.2.5b9, [example output]({{ BASE_PATH }}/assets/Pindel_version_0.2.5b9.output.zip)   
* DELLY version 0.7.6, [example output]({{ BASE_PATH }}/assets/DELLY_version_0.7.6.output.zip)  
* SVseq2 version 2.2, [example output]({{ BASE_PATH }}/assets/SVseq2_version_2.2.output.zip)  
* CNVnator version 0.3.3, [example output]({{ BASE_PATH }}/assets/CNVnator_version_0.3.3.output.zip)  
* Lumpy version 0.2.13, [example output]({{ BASE_PATH }}/assets/Lumpy_version_0.2.13.output)  
* SoftSearch version 2.4, [example output]({{ BASE_PATH }}/assets/SoftSearch_version_2.4.output)  

### intansv release version 1.7.3 support the following version of SV prediction programs
* BreakDancer version 1.4.5, [example output]({{ BASE_PATH }}/assets/BreakDancer_version_1.4.5.output)  
* Pindel version 0.2.5a8, [example output]({{ BASE_PATH }}/assets/Pindel_version_0.2.5a8.output.zip)   
* DELLY version 0.6.1, [example output]({{ BASE_PATH }}/assets/DELLY_version_0.6.1.output.zip)  
* SVseq2 version 2.2, [example output]({{ BASE_PATH }}/assets/SVseq2_version_2.2.output.zip)  
* CNVnator version 0.3, [example output]({{ BASE_PATH }}/assets/CNVnator_version_0.3.output.zip)  
* Lumpy version 0.2.8, [example output]({{ BASE_PATH }}/assets/Lumpy_version_0.2.8.output)  
* SoftSearch version 2.4, [example output]({{ BASE_PATH }}/assets/SoftSearch_version_2.4.output)  

### intansv release version 1.7.1 support the following version of SV prediction programs
* BreakDancer version 1.0r112, [example output]({{ BASE_PATH }}/assets/BreakDancer_version_1.0r112.output)  
* Pindel version 0.2.4t, [example output]({{ BASE_PATH }}/assets/Pindel_version_0.2.4t.output.zip)  
* DELLY version 0.0.11, [example output]({{ BASE_PATH }}/assets/DELLY_version_0.0.11.output.zip)  
* SVseq2 version 2.2, [example output]({{ BASE_PATH }}/assets/SVseq2_version_2.2.output.zip)  
* CNVnator version 0.2.7, [example output]({{ BASE_PATH }}/assets/CNVnator_version_0.2.7.output.zip)  
* Lumpy version 0.1.5, [example output]({{ BASE_PATH }}/assets/Lumpy_version_0.1.5.output)  
* SoftSearch version 2.3, [example output]({{ BASE_PATH }}/assets/SoftSearch_version_2.3.output)  


### Available Genome Annotation files storaged as 'Rdata' created from Gff3 files  
* [evidence-based annotation of the human genome (GRCh38), version 22 (Ensembl 79)]({{ BASE_PATH }}/assets/GffAnnotation/GRCh38.anno.RData)   
* [evidence-based annotation of the human genome (GRCh37), version 19 (Ensembl 74)]({{ BASE_PATH }}/assets/GffAnnotation/GRCh37.anno.RData)   
* [annotation of the human genome (GRCh38.p7) provided by GENCODE, Release 25]({{ BASE_PATH }}/assets/GffAnnotation/GRCh38.gencode.anno.RData)   

### Report bugs
If you find any bugs using intansv, please <a href="https://github.com/venyao/intansv/issues/intansv/issues" target="_blank">create an issue</a> or send me an email.

### Contact information
ywhzau at gmail.com


