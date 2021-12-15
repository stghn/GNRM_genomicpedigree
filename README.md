# Generalized Numerator Relationship Matrix (GNRM) for genomic and pedigree relationship matrix
## Computing Genomic Relationships according to VanRaden PM (2008)

Genomic relationship martix (GRM) was computed as:

<img src="https://render.githubusercontent.com/render/math?math=\Large G = \frac{Z Z'}{ 2 \sum pi(1-pi)}">

**Z** from the above formula was computed as: **Z = M - P**  

- where **M** is a genotype matrix of gene-content *(0/1/2)*.
   - The gene content of **0** represents when an individual for SNP *i* was major homozygous **11**
   - The gene content of **1** represents when an individual for SNP *i* was heterozygous **(12, or 21)**
   - The gene content of **2** represents when an individual for SNP *i* was minor homozygous **22**   
- where **P** is a matrix containing twice the minor allele frequency for SNP *i* (i.e. *2pi*).  

Minor allele frequencies can be computed based on  
 1.  Observed freguency from the data  
 2.  Supplied allele frequency by user 
 3.  Substituting *p=0.5* for all loci  
 4.  Average of minor allele across all loci  


## **Description of GNRM.R script**  
This R script contains the code for computing Genomic relationship matrox according to VanRaden P. (2008) and Forni et al. (2011) and computing Pedigree relationship matrix. After applying **source** function on "GNRM.R", the function called **calc_gnrm** were be available in global R environment. For computing pedigree relationship in addition to genomic relationship, you need to install and call the R-pacakge **QTLRel**.  

## **_IMPORTANT NOTE_**
 1. ONLY numeric allele codings are allowed (alleles -- 11/12/22 or genotypes -- 0/1/2).  
 2. Missing genotypes are allowed (however code missing genotypes as NA. (Adhoc imputations is done by replacing missing values with the column mean).  

**Required arguments for cal_gnrm function**
The arguments for **calc_gnrm** function are given in detail below. Please go through the information to get comprehensive explanations.

1. _genofile_: Name of the genotype file: it can be an "external file" or "R-object" (compulsory argument)
2. _genoformat_: Genotype file format. **Five** format types are allowed (**ped,tped,genotypes,Robj_ped,Robj_genotypes**) (compulsory argument)  
3. _ana_type_: **Five** options for computing GRM (**vanRaden,vanRaden_SAF,Forni_0.5,Forni_MAF,Forni_GN**) (compulsory argument)  
4. _AFREQ_: Name of a file containing all frequency (This is only used when _ana_type_ = **vanRaden_SAF**).It is an optional argument use empty string "" to represent non-use  
5. _ped_data_: Name of file with pedigree information (ID, Sire, Dam). No headers allowed. it is an optional argument use empty string "" to represent non-use  
6. _ped_option_: Specify **(TRUE)** or **(FALSE)**. when pedigree information is NOT available use F(FALSE).(compulsory argument)
7. _outputformat_: **two** output format types are allowed (**ASREML,matrix**). (compulsory argument)  
8. _outputname_: Output name of final file. This will be the text file outputted to your current directory (compulsory argument). 
9. _nIID_: Number of animals in dataset, slightly higher value increase speed of reading genotype data (compulsory argument)  
10. _missinggeno_: **TRUE** or **FALSE** if there are missing genotypes, missing genotype should be NA. (compulsory argument)  
11. _plots_: **TRUE** or **FALSE**. (compulsory argument)  

## **Some specific explanation for Arguments using cal_gnrm function**  
* Explanation for argument **_genoformat_** : This argument allows for different file format to be specified  
    1. **ped** -- PLINK ped file format, or linkage file format see PLINK.  
        format:   FamID ID sire dam sex pheno SNP1_allele1 SNP1_allele2 SNP2_allele1 SNP2_allele2 ...  

    2. **tped** -- PLINK tranposed file format, see PLINK for details. when tped format is used, specify the names of the tped file and tfam file (eg c("file.tped","file.tfam"))  
        format:   CHR SNPname geneticposition physical position sample1_SNP_allele1 sample1_SNP1_allele2 sample2_SNP_allele1 sample1_SNP_allele2 ...   

    3. **genotypes** -- SNPs are coded in genotype format (0,1,2); representing AA/AB/BB or 11/12/22  
        format:   FamID ID sire dam sex pheno SNP1 SNP2 ...  

    4. **Robj_ped** -- R-object, when you already have imported the dataset in R. Depending on the format, if it an R-object in **ped** style then use **Robj_ped** 
        format: **Robj_ped**  -- ID SNP1_allele1 SNP1_allele2 SNP2_allele1 SNP2_allele2 ...  

    5. **Robj_genotypes** -- R-object, when you already have imported the dataset in R. When the format of the R-object is in **genotypes** style then use **Robj_genotypes**.  
        format: **Robj_genotypes**  -- ID SNP1 SNP2 ...  

* Explanation for argument **_ana_type_** : Five methods in computing genomic relationships are implemented.  
    1. **vanRaden** -- allele frequency for centering and scaling genotypes are `computed from the available data for each loci`   
    2. **vanRaden_SAF** -- allele frequency for centering and scaling genotypes are `supplied by the user`. Note that you are responsible for making sure the list of allele frequencies correspond to the columns in the data file. When  **vanRaden_SAF** is used, supply allele frequency with the argument `AFREQ`.  
    3. **Forni_0.5** -- allele frequency for centering and scaling genotypes are `assumed to be 0.5 for all loci`.  
    4. **Forni_MAF** -- allele frequency for centering and scaling genotypes are `assumed to be average minor allele frequency for all loci`. The function first compute the minor allele frequency for each loci and takes the average and uses it for centering and scaling genotypes  
    5. **Forni_GN** -- allele frequency for centering and scaling genotypes are `computed from the data for each loci`. however, scaling of the genotypes are based on the trace of the `ZZ'` matrix divided by the number of genotypes individuals. read Forni et al. (2011) for more details


* Explanation for arguments **_ped_data_** and **_ped_option_**: Pedigree relationships can also be computed when the user supplies a pedigree file.  
    1. **_ped_data_** -- the user can supply a pedigree file with **three columns** only (id, sire, dam). Columns HEADERS are not allowed.  
    2. **_ped_option_** -- This argument is compulsory. The user has to supply **(TRUE)** or **(FALSE)** argument for the script to work.  


* Explanation for argument **_outputformat_**: Two output formats are allowed.  
    1. **ASREML** -- ASREML for a relationship matrix (free flow format, only the diagonal element and the lower triangle is present) - fast to export  
    2. **matrix** -- pairwise relationship in a matrix format. - faster to export  


## Sequential explanation of how to implement vanRaden (2008) in R  
 ```R
 ## Read in a genotype file (only genotypes and are coded as 0,1,2)  
 M <- read.table("example/ex_1k.genotype")[,-c(1:6)]  

## Calculate allele frequency of the second allele (i.e. minor allele frequency of genotype coded as 2)  
## Simple example for just 1 SNP with 5 genotype animals
## Example c(AA, AB, AB, BB, BB)
## Minor allele frequency of p (i.e. allele 2); sum the values in SNP1 divided by 2*nrow(SNP1)  
## This expression is equivalent to taking the mean of the column and dividing it by 2  

SNP1 <- matrix(c(0,1,1,2,2))
p <- mean(SNP1)/2
q <- 1-p   

## Calculate p and q for a geno data file using the apply function  
p <- (apply(M,2,mean))/2 
q <- 1-p  
pt2 <- 2*p  

## Subtract 2*p from M  
Z <- t(apply(M,1,function(x) x-pt2))  

## Calculate scaler K  
K <- 2*sum(p*q)  

## Compute G  
G <- (Z %*% t(Z))/K

```

## Optimise way to implement vanRaden formulae in R  
```R
# Read in a genotype file (only genotypes and are coded as 0,1,2)  (only genotypes are needed delete non-important columns) 
M <- read.table("example/ex_1k.genotype")[,-c(1:6)]
M <- scale(x=M,center=T,scale=F)
K <-sum(apply(X=M,FUN=var,MARGIN=2))  
G <- tcrossprod(M)/K  

## Examples
## Let use the script to try out some examples  
## Set working directory to the correct path  
## soucre the file  
setwd("~/packages/script_GRM/")  
source("GNRM.R")
library("QTLRel")  

## for example 1 (A) 
 # computing GRM based using PLINK - PED file format as input marker data  
 # 1k dataset (output results in _matrix_ format)  
 # ana_type is vanRaden (2008) [ vanRaden ]   

ex1mat_Gvan<-calc_gnrm(genofile="example/ex_1k.ped",genoformat="ped",ana_type="vanRaden",AFREQ="",ped_data="",ped_option=F,outputformat="matrix",outputname="ex1mat_1kGvan",nIID=300,missinggeno=F,plots=T)  

## for example 1 (B)
 # computing GRM based using PLINK - PED file format as input marker data  
 # 1k dataset (output results in _ASREML_ format)  
 # ana_type is vanRaden (2008) [ vanRaden ]   
ex1asreml_Gvan <- calc_gnrm(genofile="example/ex_1k.ped",genoformat="ped",ana_type="vanRaden",AFREQ="",ped_data="",ped_option=F,outputformat="ASREML",outputname="ex1asreml_1kGvan",nIID=300,missinggeno=F,plots=T)  

## for example 1 (C)
 # computing GRM based using PLINK - PED file format as input marker data  
 # 5k dataset (output results in _matrix_ format)
 # ana_type is vanRaden (2008) [ vanRaden ]   

ex1mat_Gvan <- calc_gnrm(genofile="example/ex_5k.ped",genoformat="ped",ana_type="vanRaden",AFREQ="",ped_data="",ped_option=F,outputformat="matrix",outputname="ex1mat_5kGvan",nIID=300,missinggeno=F,plots=T)  

## for example 1 (D)
 # Using other ana_type (grm options) 
 # Genotype are in PLINK - PED file format as input marker data 
 # 5k dataset (output results in _matrix_ format)
 # ana_type is Forni et al. (2011) [ Forni_MAF ]  

ex1_GforniMAF <- calc_gnrm(genofile="example/ex_5k.ped",genoformat="ped",ana_type="Forni_MAF",AFREQ="",ped_data="",ped_option=F,outputformat="matrix",outputname="ex1mat_5kGforniMAF",nIID=300,missinggeno=F,plots=T)  

## for example 1 (E)
 # Genotype are in PLINK - PED file format as input marker data 
 # 5k dataset (output results in _matrix_ format)
 # ana_type is Forni et al. (2011) [ Forni_GN ] 
 
ex1mat_GforniGN <- calc_gnrm(genofile="example/ex_5k.ped",genoformat="ped",ana_type="Forni_GN",AFREQ="",ped_data="",ped_option=F,outputformat="matrix",outputname="ex1mat_5kGforniGN",nIID=300,missinggeno=F,plots=T)  

## for example 2 (A)  
## computing Genomic relationship and Pedigree relationship 
 # using PLINK - PED file format as input marker data
 # using pedigree information (_ped_option=T_)
 # 5k dataset (output results in _ASREML_ format)
 # ana_type is vanRaden (2008) [ vanRaden ]
 
ex2asreml_GNRMvan <- calc_gnrm(genofile="example/ex_5k.ped",genoformat="ped",ana_type="vanRaden",AFREQ="",ped_data="example/dPedigree.txt",ped_option=T,outputformat="ASREML",outputname="ex2asreml_5kGNRMvan",nIID=300,missinggeno=F,plots=T)  

## for example 2 (B)
##computing Genomic relationship and Pedigree relationship with user supplied allele frequency
 # using PLINK - PED file format as input marker data
 # using pedigree information (_ped_option=T_)
 # 5k dataset (output results in _ASREML_ format)
 # ana_type is vanRaden (2008) [ vanRaden_SAF ]

ex2asreml_GNRMvanSAF <- calc_gnrm(genofile="example/ex_5k.ped",genoformat="ped",ana_type="vanRaden_SAF",AFREQ="example/freq_5k.txt",ped_data="example/dPedigree.txt",ped_option=T,outputformat="ASREML",outputname="ex2asreml_5kGNRMvanSAF",nIID=300,missinggeno=F,plots=T) 

## for example 3 (A)  
## computing GRM using R-object in ped format as input marker data
 # using R-object - genoformat="Robj_ped"
 # 5k dataset (output results in _natrix_ format)
 # ana_type is Forni et al. (2011) [ Forni_0.5 ]
## Steps: reading genotype file in R, and delete redundent columns except ID and specify this file for the function  

geno <- read.table("example/ex_5k.ped")[,-c(1,3:6)]
ex3mat_Gforni <- calc_gnrm(genofile="geno",genoformat="Robj_ped",ana_type="Forni_0.5",AFREQ="",ped_data="",ped_option=F,outputformat="matrix",outputname="ex3mat_5kGforni",nIID=300,missinggeno=F,plots=T)  

## for example 3 (B) 
##computing GRM using R-object in genotype format as input marker data
 # using R-object - genoformat="Robj_genotypes"
 # 1k dataset (output results in _ASREML_ format)
 # ana_type is vanRaden (2008) [ vanRaden ]
## Steps: reading genotype file in R, and delete redundent columns except ID and specify this file for the function  

ex3asreml_Gvan <- calc_gnrm(genofile="geno",genoformat="Robj_genotypes",ana_type="vanRaden",AFREQ="",ped_data="",ped_option=F,outputformat="ASREML",outputname="ex3asreml_1kGvan",nIID=300,missinggeno=F,plots=TRUE)  


## Eigen value decomposition of G can be obtained and plots made using the `eigen function` on G  
## Heatmap of G can also be created  

eigen_G <- eigen(ex1mat_Gvan)  
pca_G <- eigen_G$vectors  
pervarPCA <- round((eigen_G$values/sum(eigen_G$values))*100,3)
plot(x=pca_G[,1],y=pca_G[,2],pch=20,col="darkblue",xlab=c(paste("PCA 1 - (",pervarPCA[1],"%)",sep="")), ylab=c(paste("PCA 2 - (",pervarPCA[2],"%)",sep=""))) 
heatmap(ex1mat_Gvan,labRow=F,labCol=F,keep.dendro=T)  


## for example 4 
# Example with larger and combined multi-breed as well as crossbred population  
# GRM will be computed and PCA's will be created 
# computing might take some time  

 # marker data in **genotype format** file format as input  
 # 5k dataset (output results in _matrix_ format)
 # ana_type is vanRaden (2008) [ vanRaden ]

ex4mat_Gvan <- calc_gnrm(genofile="example/ex_diffPOP.genotypes",genoformat="genotypes",ana_type="vanRaden",AFREQ="",ped_data="",ped_option=F,outputformat="matrix",outputname="ex4mat_Gvan",nIID=4000,missinggeno=F,plots=TRUE)    
eigen_G <- eigen(ex5_Gvan)
pca_G <- eigen_G$vectors 
pervarPCA <- round((eigen_G$values/sum(eigen_G$values))*100,3)
plot(x=pca_G[,1], y=pca_G[,2], pch=20, col="darkblue", xlab=c(paste("PCA 1 - (",pervarPCA[1],"%)",sep="")), ylab=c(paste("PCA 2 - (",pervarPCA[2],"%)",sep="")))  

```

**please report bugs to soloboan@yahoo.com**

