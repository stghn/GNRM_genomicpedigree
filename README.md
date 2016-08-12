# GNRM (genomic and pedigree relationship matrix)
## Computing Genomic Relationships according to VanRaden PM (2008)

G was computed as  
G = \frac{Z Z'}{ 2 \sum pi(1-pi)}  

Z was computed as  
 Z = M - P  

where **M** genotype matrix of gene-content _(0/1/2)_.  
it was **0** if the genotype of an individual for SNP _i_ was homozygous **11**, was **1** if heterozygous **(12, or 21)**, or **2** if the genotype was homozygous **22**.  
where **P** is a matrix containing the 2 _x_ the frequency of the second allele for SNP _i_ (i.e.   _2pi_ ).  
Allele frequencies are computed as  
 1. _from the data_  
 2. _use supplied_  
 3. _0.5 for all loci_  
 4. _average frequency of the minor allele_  


## **Description of the R-function**  
Rscript (`GNRM.R`) contains the code for the computation of pairwise Genomic relationship according to VanRaden P. (2008) and Forni et al. (2011). Pedigree  relationship as also computed. When the script is sourced, the function name is `calc_gnrm`. For computing pedigree relationship in addition to genomic relationship, you need to install and call the R-pacakge `QTLRel`.  

The function argument are given in detail below. Please try and go through the information to get comprehensive explanations.

After sourcing the Rscript (GNRM.R), **GRM** can be computed with the following R-function.  
` calc_gnrm(genofile, genoformat, ana_type, AFREQ, ped_data, ped_option, outputformat, outputname, missinggeno, nIID, plots) `.

## **_important NOTE_**:  
 1. ONLY numeric allele codings are allowed (alleles -- 11/12/22 or genotypes -- 0/1/2).  
 2. missing genotypes are allowed (however code missing genotypes as NA.  
   - adhoc imputations is done by replacing missing values with the column mean).  

**R Argument**  
1. _genofile_       === (compulsory)... Name of the genotype file: it can be an "external file" or "R-object"  
2. _genoformat_     === (compulsory)... genotype file format. four (5) format types are allowed (**ped,tped,genotypes,Robj_ped,Robj_genotypes**)  
3. _ana_type_       === (compulsory)... five (5) options of computing GRM (**vanRaden,vanRaden_SAF,Forni_0.5,Forni_MAF,Forni_GN**)  
4. _AFREQ_          === (optional)... Name of a file containing all frequency (This is only used when _ana_type_ = **vanRaden_SAF**), use empty string "" to represent non-use  
5. _ped_data_       === (optional)... Name of file with pedigree information (ID, Sire, Dam). No headers allowed. it is an optional argument use empty string "" to represent non-use  
6. _ped_option_     === (compulsory)... specify **T (TRUE)** or **F (FALSE)**. when pedigree information is NOT available use F(FALSE)  
7. _outputformat_   ===  (compulsory)... three (3) output format types are allowed (**ASREML,dense,matrix**)  
8. _outputname_     === (compulsory)... output name of final file. This will be the text file outputted to your current directory  
9. _nIID_           === (compulsory)... number of animals in dataset, slightly higher value increase speed of reading genotype data  
10. _missinggeno_   === (compulsory)... **TRUE/T** or **FALSE/F** if there are missing genotypes, missing genotype should be NA  
11. _plots_          === (compulsory)...  **TRUE/T** or **FALSE/F**  

## **Some specific explanation of Argument 2, 3, 4, 5 and 7**  
* Explanation for argument 2 [_genoformat_] : genoformat allows for different file format to be specified  
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

* Explanation for argument 3 [_ana_type_] : Five (5) methods in computing genomic relationships are implemented.  
These methods are based on vanRaden (2008) and Forni et al. (2011). The main differences between the methods are which allele frequencies are used.  

    1. **vanRaden** -- allele frequency for centering and scaling genotypes are `computed from the available data for each loci`   
    2. **vanRaden_SAF** -- allele frequency for centering and scaling genotypes are `supplied by the user`. Note that you are responsible for making sure the list of allele frequencies correspond to the columns in the data file. When  **vanRaden_SAF** is used, supply allele frequency with the argument `AFREQ`.  
    3. **Forni_0.5** -- allele frequency for centering and scaling genotypes are `assumed to be 0.5 for all loci`.  
    4. **Forni_MAF** -- allele frequency for centering and scaling genotypes are `assumed to be average minor allele frequency for all loci`. The function first compute the minor allele frequency for each loci and takes the average and uses it for centering and scaling genotypes  
    5. **Forni_GN** -- allele frequency for centering and scaling genotypes are `computed from the data for each loci`. however, scaling of the genotypes are based on the trace of the `ZZ'` matrix divided by the number of genotypes individuals. read Forni et al. (2011) for more details


* Explanation for argument 4 and 5 [_ped_data_ , _ped_option_] : Pedigree relationships can also be computed when the user supplies a pedigree file.  
    **_ped_data_** -- the user can supply a  pedigree file with three (3) columns only (id, sire, dam). Columns HEADERS are not allowed.  
    **_ped_option_** -- This argument is compulsory. The user has to supply a T (TRUE) or F (FALSE) argument for the script to work.  


* Explanation for argument 7 [_outputformat_] : Three (3) output format are allowed.  
    1. **ASREML** -- ASREML for a relationship matrix (free flow format, only the diagonal element and the lower triangle is present) - fast to export  
    2. **dense**  -- free flow format with all pairwise relationship present - slow to export  
    3. **matrix** -- pairwise relationship in a matrix format. - faster to export  


## Sequential explanation of how to implement vanRaden (2008) in R  
read in a genotype file (only genotypes and are coded as 0,1,2)  
`M <- read.table("example/ex_1k.genotype")[,-c(1:6)]`  

calculate allele frequency of the second allele (i.e. allele frequency of genotype coded as 2)  
simple example for just 1 SNP with 5 genotype animals
Example c(AA, AB, AB, BB, BB)
allele frequency of p (i.e. allele 2); sum the values in SNP1 divided by 2*nrow(SNP1)  
This expression is equivalent to taking the mean of the column and dividing it by 2  

`SNP1 <- matrix(c(0,1,1,2,2))`  
`p <- mean(SNP1)/2`   
`q <- 1-p`   

calculate 'p' and 'q' for a geno data file using the apply function  
`p <- (apply(M,2,mean))/2`  
`q <- 1-p`  
`pt2 <- 2*p`  

###### subtract 2*p from M  
Z <- t(apply(M,1,function(x) x-pt2))  

###### calculate scaler K  
K <- 2*sum(p*q)  

###### compute G  
G <- (Z %*% t(Z))/K  


##### Optimise way to implement vanRaden formulae in R  

##### read in a genotype file (only genotypes and are coded as 0,1,2)  
M <- read.table("example/ex_1k.genotype")[,-c(1:6)]  
###### first six non-important columns deleted  
M <- scale(x=M,center=T,scale=F)  
K<-sum(apply(X=M,FUN=var,MARGIN=2))  
G <- tcrossprod(M)/K  


#### Let use the script to try out some examples  
##### Set working directory to the correct path  
setwd("~/packages/script_GRM/")  

#### soucre the file  
source("GNRM.R")  
library("QTLRel")  

** _for example 1_**  
computing GRM based using PLINK - PED file format as input marker data  
output GRM format is '_matrix_' type and '_ASREML_'  


#### 1k dataset (output results in matrix format)  
##### ana_type is vanRaden (2008)  
ex1mat_Gvan <- calc_gnrm(genofile="example/ex_1k.ped",genoformat="ped",ana_type="vanRaden",AFREQ="",ped_data="",ped_option=F,outputformat="matrix",outputname="ex1mat_1kGvan",nIID=300,missinggeno=F,plots=T)  

#### 1k dataset (output results in ASREML format)  
##### ana_type is vanRaden (2008)  
ex1asreml_Gvan <- calc_gnrm(genofile="example/ex_1k.ped",genoformat="ped",ana_type="vanRaden",AFREQ="",ped_data="",ped_option=F,outputformat="ASREML",outputname="ex1asreml_1kGvan",nIID=300,missinggeno=F,plots=T)  

#### Try the 5k dataset (output results in matrix format)  
##### ana_type is vanRaden (2008) [ vanRaden ]   
ex1mat_Gvan <- calc_gnrm(genofile="example/ex_5k.ped",genoformat="ped",ana_type="vanRaden",AFREQ="",ped_data="",ped_option=F,outputformat="matrix",outputname="ex1mat_5kGvan",nIID=300,missinggeno=F,plots=T)  


#### Example for using other ana_type (grm options)  
##### ana_type is Forni et al. (2011) [ Forni_MAF ]  
##### output as Matrix format  
ex1_GforniMAF <- calc_gnrm(genofile="example/ex_5k.ped",genoformat="ped",ana_type="Forni_MAF",AFREQ="",ped_data="",ped_option=F,outputformat="matrix",outputname="ex1mat_5kGforniMAF",nIID=300,missinggeno=F,plots=T)  

  # ana_type is Forni et al. (2011) [ Forni_GN ]  
  # output as Matrix format  
ex1mat_GforniGN <- calc_gnrm(genofile="example/ex_5k.ped",genoformat="ped",ana_type="Forni_GN",AFREQ="",ped_data="",ped_option=F,outputformat="matrix",outputname="ex1mat_5kGforniGN",nIID=300,missinggeno=F,plots=T)  


**_for example 2_**:  
computing Genomic relationship and Pedigree relationship using PLINK - PED file format as input marker data output GRM formats _ASREML_  

###### based on vanRaden (2008) with allele frequency computed from the data 
ex2asreml_GNRMvan <- calc_gnrm(genofile="example/ex_5k.ped",genoformat="ped",ana_type="vanRaden",AFREQ="",ped_data="example/dPedigree.txt",ped_option=T,outputformat="ASREML",outputname="ex2asreml_5kGNRMvan",nIID=300,missinggeno=F,plots=T)  

###### based on vanRaden (2008) but with user supplied allele frequency for each loci [ vanRaden_SAF ]  
ex2asreml_GNRMvanSAF <- calc_gnrm(genofile="example/ex_5k.ped",genoformat="ped",ana_type="vanRaden_SAF",AFREQ="example/freq_5k.txt",ped_data="example/dPedigree.txt",ped_option=T,outputformat="ASREML",outputname="ex2asreml_5kGNRMvanSAF",nIID=300,missinggeno=F,plots=T)  


**_for example 3_**:  
computing GRM using genotype file format as input marker data output GRM format is '_matrix_' type  

##### ana_type is vanRaden (2008) [ vanRaden ] 
##### output as Matrix format  
ex3mat_Gvan <- calc_gnrm(genofile="example/ex_1k.genotype",genoformat="genotypes",ana_type="vanRaden",AFREQ="",ped_data="",ped_option=F,outputformat="matrix",outputname="ex3asreml_1kGvan",nIID=300,missinggeno=F,plots=T)  

##### ana_type is Forni et al. (2011) [ Forni_0.5 ]  
##### output as Matrix format  
ex3mat_Gforni <- calc_gnrm(genofile="example/ex_1k.genotype",genoformat="genotypes",ana_type="Forni_0.5",AFREQ="",ped_data="",ped_option=F,outputformat="matrix",outputname="ex3mat_1kGforni",nIID=300,missinggeno=F,plots=T)  


**_for example 4_**:  
computing GRM using **R-object in ped format** as input marker data or using **R-object in genotype format** as input marker data output GRM format is '_matrix_' type  

##### R-object in ped format  
###### reading genotype file in R,  
###### delete redundent columns except ID and specify this file for function  
geno <- read.table("example/ex_5k.ped")[,-c(1,3:6)]  

##### ana_type is Forni et al. (2011) [ Forni_0.5 ] 
##### output as Matrix format  
ex4mat_Gforni <- calc_gnrm(genofile="geno",genoformat="Robj_ped",ana_type="Forni_0.5",AFREQ="",ped_data="",ped_option=F,outputformat="matrix",outputname="ex4mat_5kGforni",nIID=300,missinggeno=F,plots=T)  


##### **R-object genotype format**  
##### reading genotype file in R, 
##### delete redundent columns except ID and specify this file for function  
geno <- read.table("example/ex_1k.genotype")[,-c(1,3:6)]  

##### ana_type is vanRaden (2008) [ vanRaden ] 
##### output as ASREML format  
ex4asreml_Gvan <- calc_gnrm(genofile="geno",genoformat="Robj_genotypes",ana_type="vanRaden",AFREQ="",ped_data="",ped_option=F,outputformat="ASREML",outputname="ex4asreml_1kGvan",nIID=300,missinggeno=F,plots=TRUE)  


Eigen value decomposition of G can be obtained and plots made using the `eigen function` on **G**  
Heatmap of G can also be created  

eigen_G <- eigen(ex1mat_Gvan)  
pca_G <- eigen_G$vectors  
pervarPCA <- round((eigen_G$values/sum(eigen_G$values))*100,3)  
plot(x=pca_G[,1],y=pca_G[,2],pch=20,col="darkblue",xlab=c(paste("PCA 1 - (",pervarPCA[1],"%)",sep="")), ylab=c(paste("PCA 2 - (",pervarPCA[2],"%)",sep="")))  
heatmap(ex1mat_Gvan,labRow=F,labCol=F,keep.dendro=T)  


**Example with larger and combined multi-breed as well as crossbred population**
_GRM will be computed and PCA's will be created_ 
computing might take some time  
computing GRM using **genotype format** as input marker data  
output GRM format is '_matrix_' type  

ex5mat_Gvan <- calc_gnrm(genofile="example/ex_diffPOP.genotypes",genoformat="genotypes",ana_type="vanRaden",AFREQ="",ped_data="",ped_option=F,outputformat="matrix",outputname="ex5mat_Gvan",nIID=4000,missinggeno=F,plots=TRUE)  
eigen_G <- eigen(ex5_Gvan)  
pca_G <- eigen_G$vectors  
pervarPCA <- round((eigen_G$values/sum(eigen_G$values))*100,3)  
plot(x=pca_G[,1], y=pca_G[,2], pch=20, col="darkblue", xlab=c(paste("PCA 1 - (",pervarPCA[1],"%)",sep="")), ylab=c(paste("PCA 2 - (",pervarPCA[2],"%)",sep="")))  

**please report bugs to soloboan@yahoo.com**

