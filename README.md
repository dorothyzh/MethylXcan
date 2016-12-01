README for MethylXcan
===============

[Huan ZHONG](https://github.com/dorothyzh/) \(hzhong5 at uab dot edu\)

* * *

Table of Contents
-----------------
* [Introduction](#introduction)
* [Prerequisites & Installation](#compilation)
* [Usage](#usage)
    * [Using Methylation to predict Gene expression](#built)
* [Demo](#demo)
* [Authors](#authors)

* * *

## <a name="introduction"></a> Introduction
We have developed the MethylXcan which can predict the gene expression pattern based on DNA methylation data, and can help to provide insights into the mechanism of these associations.

## <a name="compilation"></a> Prerequisites & Installation
R 3.2.1 is suggested. Some R packages, like "glmnet" and "methods" are also required.They will be automatically installed.
No further installation is needed. You only need to format the input files acording to the requirement, and run one perl script on these files.

## <a name="usage"></a> Usage

### I. Preparing Input Files

#### a) ex_probe_list.txt
   One tab-delimited annotation file containing gene expression probe, gene name, official name, chromosome and locations. Here is one gene entry as example.

        
         ILMN_2038774    EEF1A1          NM_001402.5     chr6:74284964-74285013
       
    **Note**  This file is suggested without header.

#### b) ex_dataset.txt 
   One tab-delimited gene expression profiling dataset, containing gene probe and its profiling values (normalized if it is microarray data) from different samples.
       
       Hybridization REF       TWPID6598       TWPID3283       TWPID5553...
       ILMN_1343291    16.043236443862 15.9458304153505        15.9085900238073...
       
   
       

#### c) me_dataset.txt 
   One tab-delimited DNA methylation dataset, containing CpG probes and their methylation values (normalized if it is microarray data, beta values are required) from different samples.
       
       Hybridization REF       TWPID5259       TWPID8404       TWPID2116...
       cg00240178      0.36676 0.38544 0.30756...
       
   
       

#### d) methylation_annotation.txt
   One tab-delimited CpG probe annotation file.
   
       IlmnID  CHR     MAPINFO Strand  UCSC_RefGene_Name       UCSC_RefGene_Group
       cg00240178      6       74232108        R       EEF1A1  TSS1500
       
#### e) gene_annotation.txt
   One tab-delimited gene annotation file.
   
        chr     strand  txStart txEnd   name
        chr6    -       74225472        74230755        EEF1A1


### II. Run the perl script on the prepared input files.
      
        perl run_gene_list.pl  ex_probe_list.txt  ex_dataset.txt  me_dataset.txt methylation_annotation.txt  data/gene_annotation.demo.txt

### III.  Final results.
   The final results will be named as "MethylXcan.txt", including 21 columns.


__CpG:__ name of CpG probes.  

__n.site:__ number of CpG sites per gene.


__gene:__ gene name.

__beta.single:__  regression coefficient from single regression of gene expression on its each CpGs methylation separately.

__beta.multiple:__ regression coefficients from multiple regression of gene expression on the methylation of its all CpG sites simultaneously. 

__beta.glmnet:__ coefficient from lasso regression between gene expression and its corresponding CpGs' methylation ratios.

__R2.single.max:__  the largest coefficient of determination from the single regressions of one gene.

__R2.single.var:__  the variance of coefficient of determination from the single regressions of one gene.

__R2.single.cv.max:__ max coefficient of determination from cross-validation of single regression.

__R2.single.cv.max.var:__ variance between coefficients of determination from cross-validation of single regression.

__R2.multiple:__ coefficient of determination from multiple regressions.

__R2.multiple.adjust:__ adjusted coefficient of determination from multiple regressions.

__R2.multiple.cv:__ coefficient of determination from cross-valudation of multiple regressions.

__R2.multiple.cv.var:__ variance of coefficient of determination from cross-valudation of multiple regressions.

__R2.glmnet:__ coefficient of determination from lasso regressions.

__R2.glmnet.cv:__ coefficient of determination from cross-validation of lasso regressions.

__R2.glmnet.cv.var:__ variance of coefficient of determination from cross-validation of lasso regressions.

__p.single:__ p-value from single regression.

__p.multiple:__ p-value for each CpG in a multiple regressions.

__p.multiple.overall:__ the overall p-value from multiple regressions.

__genevar:__ variance of gene expression profiling between different samples.

__dist:__ the distance between each CpG and its corresponding gene's TSS site.


### IV. Computing Time.
The program might take a long time to run, hours for Gb-sized datasets. In demo, it might take 10 seconds to run 4 probes. So when running the job in cluster, it is recommended to split your probe files (ex_probe_list.txt) into several files, and send the jobs to different nodes. 



## <a name="demo"></a> Demo

Download the demo folder, and go into the demo folder and simply run 
   
    perl script/run_gene_list.pl \
            data/ex_probe_list.demo.txt \
            data/ex_dataset.demo.txt \
            data/me_dataset.demo.txt \
            data/methylation_annotation.demo.txt \
            data/gene_annotation.demo.txt

The final "MethylXcan.txt" is the final results.


            
## <a name="authors"></a> Authors
............

