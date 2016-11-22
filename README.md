README for MethylXcan
===============

[Huan ZHONG](https://github.com/dorothyzh/) \(zhdorothy5 at uab dot edu\)

* * *

Table of Contents
-----------------
* [Introduction](#introduction)
* [Compilation & Installation](#compilation)
* [Usage](#usage)
    * [Using Methylation to predict Gene expression](#built)
* [Demo](#demo)
* [Authors](#authors)

* * *

## <a name="introduction"></a> Introduction
We have developed the MethylXcan which can predict the transcriptome profiling lanscape based on DNA methylation data, and provide insights into the mechanism of these associations.

## <a name="compilation"></a> Prerequisites & Installation
R 3.2.1 is suggested. Some R packages, like "glmnet" and "methods" are also required.
No installation should be performed. You only need to prepare the proper format of input files, and run one perl script on these files.

## <a name="usage"></a> Usage

### I. Preparing Input Files

#### a) ex_probe_list.txt
   One tab delimited annotation file, containing gene expression probe, gene name, official name, chromosome and locations.

       ILMN_2038774    EEF1A1  NM_001402.5     chr6:74284964-74285013
       
    **Note** that this file could contain no header.

#### b) ex_dataset.txt 
   One tab delimited gene expression profiling dataset, containing gene probe and its profiling values(normalized if it is microarray data) from different samples.
       
       Hybridization REF       TWPID6598       TWPID3283       TWPID5553...
       ILMN_1343291    16.043236443862 15.9458304153505        15.9085900238073...
       
   
       

#### c) me_dataset.txt 
   One tab delimited DNA methylation dataset, containing cpg probe and its methylation values(normalized if it is microarray data, beta values are required) from different samples.
       
       Hybridization REF       TWPID5259       TWPID8404       TWPID2116...
       cg00240178      0.36676 0.38544 0.30756...
       
   
       

#### d) methylation_annotation.txt
   One tab delimited cpg probe annotation
   
       IlmnID  CHR     MAPINFO Strand  UCSC_RefGene_Name       UCSC_RefGene_Group
       cg00240178      6       74232108        R       EEF1A1  TSS1500
       
#### e) gene_annotation.txt
   One tab delimited gene annotation.
   
        chr     strand  txStart txEnd   name
        chr6    -       74225472        74230755        EEF1A1


### II. Run the perl script on the prepared input files.
      
        perl run_gene_list.pl  ex_probe_list.txt  ex_dataset.txt  me_dataset.txt methylation_annotation.txt  data/gene_annotation.demo.txt

#### III. Check the final results.
   The final results will be named as "MethylXcan.txt", including 22 columns.
   colname = c("CpG","n.site","gene","beta.single","beta.multiple","beta.glmnet",
              "R2.single.max","R2.single.max.var","R2.single.cv.max","R2.single.cv.max.var",
              "R2.multiple","R2.multiple.adjust","R2.multiple.cv","R2.multiple.cv.var",
              "R2.glmnet","R2.glmnet.cv","R2.glmnet.cv.var",
              "p.single","p.multiple","p.multiple.overall","genevar","dist"
 

## <a name="demo"></a> demo

Download the demo folder, and go into the demo folder and simply run 
   
    perl script/run_gene_list.pl \
            data/ex_probe_list.demo.txt \
            data/ex_dataset.demo.txt \
            data/me_dataset.demo.txt \
            data/methylation_annotation.demo.txt \
            data/gene_annotation.demo.txt

The final "MethylXcan.txt" is the final results.


            



