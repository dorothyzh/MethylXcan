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

#### A.ex_probe_list.txt
A tab delimited annotation file, containing probe name, gene name, official name, chromosome and locations.

      ILMN_2038774    EEF1A1  NM_001402.5     chr6:74284964-74285013

colname = c("CpG","n.site","gene","beta.single","beta.multiple","beta.glmnet",
              "R2.single.max","R2.single.max.var","R2.single.cv.max","R2.single.cv.max.var",
              "R2.multiple","R2.multiple.adjust","R2.multiple.cv","R2.multiple.cv.var",
              "R2.glmnet","R2.glmnet.cv","R2.glmnet.cv.var",
              "p.single","p.multiple","p.multiple.overall","genevar","dist"


## <a name="demo"></a> demo

Download the demo folder, and go into the demo folder and simply run 
   
    perl script/run_gene_list.pl \
            data/ex_probe_list.demo.txt \
            data/me_probe_annotation.demo.csv \
            data/ex_dataset.demo.txt \
            data/me_dataset.demo.txt \
            data/methylation_annotation.demo.txt \
            
            
            
            
            
            data/gene_annotation.demo.txt

The final "MethylXcan.txt" is the final results.


            



