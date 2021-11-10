# Welcome to CropGBM！

[![DOI](https://zenodo.org/badge/243508749.svg)](https://zenodo.org/badge/latestdoi/243508749) [![Anaconda-Server Badge](https://anaconda.org/xu_cau_cab/cropgbm/badges/platforms.svg)](https://anaconda.org/xu_cau_cab/cropgbm) [![Anaconda-Server Badge](https://anaconda.org/xu_cau_cab/cropgbm/badges/installer/conda.svg)](https://conda.anaconda.org/xu_cau_cab)

## Introduction

Crop Genomic Breeding machine (CropGBM) is a multifunctional Python3 program that integrates data preprocessing, population structure analysis, SNP selection, phenotype prediction, and data visualization. Has the following advantages:

* Use LightGBM algorithm to quickly and accurately predict phenotype values and support GPU-accelerated training.
* Supports selection and visualization of SNPs that are strongly related to phenotype.
* Support PCA and t-SNE two dimensionality reduction algorithms to extract SNP information.
* Support Kmeans and OPTICS two clustering algorithms to analyze the sample population structure.
* Plot histograms of heterozygosity rate, deletion rate, and frequency of alleles for genotype data.


## Documentation

*English version documentation*: [https://ibreeding.github.io](https://ibreeding.github.io)

*Chinese version documentation*: [https://ibreeding-ch.github.io](https://ibreeding-ch.github.io)


## Download

Download source code : [https://github.com/YuetongXU/CropGBM/releases/tag/cropgbm-v1.1.2](https://github.com/YuetongXU/CropGBM/releases/tag/cropgbm-v1.1.2)


## Installation

### Install via Conda (Recommend)

    $ conda install -c xu_cau_cab cropgbm 

### Install via pip

    $ pip install --user cropgbm

### Install via source code

    $ tar -zxf CropGBM.tar.gz

    # Install Python package dependencies of CropGBM: setuptools, wheel, numpy, scipy, pandas, scikit-learn, lightgbm, matplotlib, seaborn
    $ pip install --user setuptools wheel numpy scipy pandas scikit-learn lightgbm matplotlib seaborn
    
    # Install external dependencies of CropGBM: PLINK 1.90 
    $ wget s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20191028.zip
    $ mkdir plink_1.90
    $ unzip plink_linux_x86_64_20191028.zip -d ./plink_1.90
    
    # Add CropGBM, PLINK to the system environment variables for quick use:
    $ vi ~/.bashrc
    export PATH="/userpath/CropGBM:$PATH"
    export PATH="/userpath/plink1.90:$PATH"
    $ source ~/.bashrc


## Test (For Conda)

Enter the ‘/miniconda3/pkgs/cropgbm-1.1.2-py39_0/info/test’  folder
Run the `run_test.py` to check whether cropgbm can run successfully locally.



## About

**Citation**: Jun Yan, Yuetong Xu, Qian Cheng, Shuqin Jiang, Qian Wang, Yingjie Xiao, Chuang Ma, Jianbing Yan and Xiangfeng Wang. _LightGBM: accelerated genomically-designed crop breeding through ensemble learning._ 

**Supplementary Information**: Support data and materials for the manuscript is available at _https://github.com/YuetongXU/Cropgbm-Paper_

**Contact us**: cropgbm@163.com

**Note**: Academic users can download directly, industrial users first contact us.





