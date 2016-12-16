calc_gnrm <- function(genofile,genoformat,ana_type,AFREQ,ped_data,ped_option,outputformat,outputname,missinggeno,nIID,plots){
  ##### ATTENTION ....... PLEASE SCROLL DOWN TO READ ABOUT THE SCRIPT   ######
  
  ######################################################################################
  if(genoformat!="tped"){
    if(missing(genofile)){
      stop(cat("\n  Need to specify the input genotype file , it can be an \n 'external file' or 'R-object'"))
    }
  } else if(genoformat=="tped"){
    genofile1 <- genofile[1]
    genofile2 <- genofile[2]
    if(missing(genofile1)){
      stop(cat("\n  Need to specify the input genotype file , it can be an \n 'external file' or 'R-object'"))
    }
    if(missing(genofile2)){
      stop(cat("\n  Need to specify the input genotype file , it can be an \n 'external file' or 'R-object'"))  
  }
  }
  
  if(missing(genoformat))
    stop(cat("\n  Need to specify genotype file format ..... (choose one of them)\n'1. ped'\n'2. tped'\n'3. genotypes'\n'4. Robj_ped'\n'5. Robj_genotypes'"))
  if(missing(ana_type))
    stop('\n  Specify the grm option you want to use  ..... (choose one of them) \n"1. vanRaden"\n"2. vanRaden_SAF"\n"3. Forni_0.5"\n"4. Forni_MAF"\n"5. Forni_GN"')
  if(missing(outputformat))
    stop("\n  Specify the format of the output file ..... (choose one of them)\n'1. ASREML'\n'2. dense'\n'3. matrix' ")
  if(missing(outputname))
    stop("\n  Specify the name of the output file ")
  if(missing(nIID))
    stop("\n Specify the minimum number of animals in the dataset")
  if(missing(missinggeno))
    stop("\n Specify if missing data is in dataset")
  
  #############################################################################
  gF <- c("ped","tped","genotypes","Robj_ped","Robj_genotypes")
  oF <- c("ASREML","dense","matrix")
  opt <- c("vanRaden","vanRaden_SAF","Forni_0.5","Forni_MAF","Forni_GN")
  
  if(!genoformat %in% gF)
    stop(cat("\n Specify the genoformat correctly\n'1. ped'\n'2. tped'\n'3. genotypes'\n'4. Robj_ped'\n'5. Robj_genotypes'"))
  if (!outputformat %in% oF)
    stop(cat("\n Specify the outputformat correctly \n'1. ASREML'\n'2. dense'\n'3. matrix'"))
  if(!ana_type %in% opt)
    stop(cat('\n Specify the correct ana_type option \n one of the following ...',
             '\n"1. vanRaden"\n"2. vanRaden_SAF"\n"3. Forni_0.5"\n"4. Forni_MAF"\n"5. Forni_GN",\n'))
  
  ##########################################################################################
  
  #checking if files specified exist
  if(genoformat=="ped" | genoformat=="genotypes"){
    if(file.exists(genofile)==F)
      stop(cat(paste(genofile," does not exit, re-specify the correct path or filename",sep="")))
  } else if(genoformat=="Robj_ped" | genoformat=="Robj_genotypes"){
    if(exists(genofile)==F)
      stop(cat(paste(genofile," does not exist as an R-object",sep="")))
  } else if(genoformat=="tped"){
    if(file.exists(genofile[1])==F)
      stop(cat(paste(genofile[1]," does not exit, re-specify the correct path or filename",sep="")))
    if(file.exists(genofile[2])==F)
      stop(cat(paste(genofile[2]," does not exit, re-specify the correct path or filename",sep="")))
  }
  ##############################################################################
  
  # checking if file exist for computing vanRaden_SAF
  if(ana_type=="vanRaden_SAF"){
    if(file.exists(AFREQ)==F)
      stop(cat(paste(AFREQ," does not exit, re-specify the correct path or filename",sep="")))
  }
  
  #if(ped_option==T){library(newm)}
  
  #############################################################################
  if (genoformat=="ped"){
    cat('\n....... Importing genotype file ..........\n')
    nc <- ncol(read.table(genofile,header=F,nrows=2)) 
    geno  <- read.table(genofile,header=F,nrows=nIID,colClasses=c(rep("character",6),rep("numeric",nc-6)))
    cat('....... Genotype file imported ..........\n\n')
    cat('....... file processing for GRM computation started ..........\n')
    IID   <- as.vector(geno[,2])
    geno  <- geno[,-1:-6]
    geno  <- geno-1
    geno  <- geno[,seq(1,ncol(geno),2)] + geno[,seq(2,ncol(geno),2)]  
    if(missinggeno==T){
      m.geno <-round(colMeans(geno,na.rm=TRUE),1) 
      missgeno<-function(x){z <- geno[,x]; z[is.na(z)]=m.geno[x]; return(z)}
      geno[] <-lapply(1:(ncol(geno)),missgeno)
    }
    
  } else if (genoformat=="genotypes"){
    cat('\n....... Importing genotype file ..........\n')
    nc <- ncol(read.table(genofile,header=F,nrows=2)) 
    geno  <- read.table(genofile,header=F,nrows=nIID,colClasses=c("character",rep("numeric",nc-1)))
    cat('....... Genotype file imported ..........\n\n')
    cat('....... file processing for GRM computation started ..........\n')
    IID  <- as.vector(geno[,2])
    geno <- geno[,-1:-6]
    if(missinggeno==T){
      m.geno <-round(colMeans(geno,na.rm=TRUE),1) 
      missgeno<-function(x){z <- geno[,x]; z[is.na(z)]=m.geno[x]; return(z)}
      geno[] <-lapply(1:(ncol(geno)),missgeno)
    }
    
  } else if (genoformat=="tped"){
    cat('\n....... Importing genotype file ..........\n')
    nc <- ncol(read.table(genofile[1],header=F,nrows=2))
    geno  <- read.table(paste(genofile[1],sep=""),header=F,colClasses=c(rep("character",4),rep("numeric",nc-4)))
    fam  <- read.table(paste(genofile[2],sep=""),header=F,colClasses=c(rep("character",6)))
    cat('....... Genotype file imported ..........\n\n')
    cat('....... file processing for GRM computation started ..........\n')
    IID   <- as.vector(fam[,2])
    geno  <- geno[,-1:-4]
    geno  <- geno-1
    geno  <- geno[,seq(1,ncol(geno),2)] + geno[,seq(2,ncol(geno),2)]
    geno <- t(geno)
    if(missinggeno==T){
      m.geno <-round(colMeans(geno,na.rm=TRUE),1)
      missgeno<-function(x){z <- geno[,x]; z[is.na(z)]=m.geno[x]; return(z)}
      geno[] <-lapply(1:(ncol(geno)),missgeno)
    }
    
  } else if (genoformat=="Robj_ped"){
    geno <- get(genofile)
    cat('\n....... file processing for GRM computation started ..........\n\n')
    IID <- as.vector(geno[,1])
    geno <- geno[,-1]
    geno <- geno-1
    geno <- geno[,seq(1,ncol(geno),2)] + geno[,seq(2,ncol(geno),2)]
    if(missinggeno==T){
      m.geno <-round(colMeans(geno,na.rm=TRUE),1) 
      missgeno<-function(x){z <- geno[,x]; z[is.na(z)]=m.geno[x]; z}
      geno[] <-lapply(1:(ncol(geno)),missgeno)
    }
    
  } else if (genoformat=="Robj_genotypes"){
    cat('\n....... file processing for GRM computation started ..........\n\n')
    geno <- get(genofile)
    IID <- as.vector(geno[,1])
    geno <- geno[,-1]
    if(missinggeno==T){
      m.geno <-round(colMeans(geno,na.rm=TRUE),1) 
      missgeno<-function(x){z <- geno[,x]; z[is.na(z)]=m.geno[x]; return(z)}
      geno[] <-lapply(1:(ncol(geno)),missgeno)
    }
  }
  ################################################################################
  
  # Methods used in computing GRM (basd on vanRaden (2008) and Forni et al. (2011))
  if (ana_type=="vanRaden"){
    Z <- scale(x=geno,center=T,scale=F)
    K <- sum(apply(X=geno,FUN=function(x){p=mean(x)/2;p=2*p*(1-p)},MARGIN=2))
    G <- tcrossprod(Z)/K
    cat('....... G computed according to vanRaden (2008) ..........\n')
    cat('....... Allele frequency were the observed frequency ........\n\n')
  }
  if(ana_type=="vanRaden_SAF"){
    afreq <- as.vector(read.table(AFREQ,header=F,colClasses="numeric")[,1])
    
    if(length(afreq)<ncol(geno))
      stop("number of rows of Allele frequency files less than number of SNP/genotypes")
    
    Z <- scale(x=geno,center=(2*afreq),scale=F)
    K <- sum(2*afreq*(1-afreq))
    G <- tcrossprod(Z)/K
    cat('....... G computed according to vanRaden (2008) ........ \n')
    cat('....... User supplied allele frequency  ..........\n\n')
  }
  if (ana_type=="Forni_0.5"){
    afreq <- 0.5
    Z <- scale(x=geno,center=rep(2*afreq,ncol(geno)),scale=F)
    K <- (2*afreq*(1-afreq))*ncol(geno)  
    G <- tcrossprod(Z)/K
    cat('....... G computed according to Forni (2011) ........ \n')
    cat('....... Allele frequency was fixed at 0.5 for all loci ........\n\n')
  }
  if(ana_type=="Forni_MAF"){
    afreq <- apply(X=geno,FUN=function(x){p=mean(x)/2},MARGIN=2)
    afreq <- ifelse(afreq>0.5,(1-afreq),afreq); afreq<- mean(afreq)
    Z <- scale(x=geno,center=rep(2*afreq,ncol(geno)),scale=F)
    K <-(2*afreq*(1-afreq))*ncol(geno)  
    G <- tcrossprod(Z)/K
    cat('....... G computed according to Forni (2011) ........ \n') 
    cat('....... Allele frequency fixed as the ........\n')
    cat('....... average frequency of the minor allele (MAF) .....\n') 
    cat('....... Average MAF was ---',round(afreq,4),' ........\n\n')  
  }
  if(ana_type=="Forni_GN"){
    Z <- scale(x=geno,center=T,scale=F)
    G <- tcrossprod(Z)
    K <-sum(diag(G))/nrow(geno)
    G <- G/K
    cat('....... G computed according to Forni (2011) ........ \n') 
    cat('....... Normalised matrix .......\n')
    cat('....... allele frequency was the observed frequency for each loci ......\n\n')
  }
  ######################################################################################
  
  # NRM calculation if pedigree is available
  # This might take time and computer memory
  if(ped_option==T){
    ped_data <- read.table(ped_data,header=F,colClasses=c(rep("character",3)))
    colnames(ped_data) <- c("id","sire","dam")
    fullA <- 2*(kinship(ped_data,ids=IID)) 
  }
  ########################################################################################
  
  #remove genofile and pedigree data
  rm(geno,Z,ped_data)
  
  ########################################################################################
  
  if(length(IID)>20){
    itercheck <- round(length(IID)/20,digits=0)
  } else {
    itercheck <- round(length(IID)/4,digits=0)
  }
  
  ##########################################################################################
  
  if(outputformat=="ASREML" & ped_option==F){
    cat('....... Output file preparation started ..........\n')
      Glist <- as.data.frame(which(row(G)>=col(G),arr.ind=TRUE))
      Glist$G <- G[lower.tri(G,diag=T)]
      Glist <- Glist[,c(2,1,3)]
      Glist <- Glist[order(Glist[,2],Glist[,1]),]
      Glist <- Glist[,c(2,1,3)]
    write.table(Glist,paste(outputname,".grm",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
    cat(paste('....... Output file exported ',outputname,'.grm',' ..........\n',sep=''))
  } else if(outputformat=="ASREML" & ped_option==T){
    cat('....... Output file preparation started ..........\n')
      Glist <- as.data.frame(which(row(G)>=col(G),arr.ind=TRUE))
      Glist$G <- G[lower.tri(G,diag=T)]
      Glist$A <- fullA[lower.tri(fullA,diag=T)]
      Glist <- Glist[,c(2,1,3,4)]
      Glist <- Glist[order(Glist[,2],Glist[,1]),]
      Glist <- Glist[,c(2,1,3,4)]
    write.table(Glist,paste(outputname,".grm",sep=""),col.names=F,row.names=F,quote=F,sep="\t")
    cat(paste('....... Output file exported ',outputname,'.grm',' ..........\n',sep=''))
    cat('+++++++++ output file columns ++++++++++\n')
    cat(' originalID1  originalID2 recodedID1 recodedID2 G_values Ped_values \n\n')
  } 
  
  if(outputformat=="matrix" & ped_option==F){
    cat('....... Output file preparation started ..........\n')
    colnames(G) <- IID
    rownames(G) <- IID
    write.table(G,paste(outputname,".grm",sep=""),col.names=T,row.names=F,quote=T,sep="\t")
    cat(paste('....... Output file exported ',outputname,'.grm',' ..........\n',sep=''))
  } else if(outputformat=="matrix" & ped_option==T){
    cat('....... Output file preparation started ..........\n')
    colnames(G) <- IID
    rownames(G) <- IID
    write.table(G,paste(outputname,".grm",sep=""),col.names=T,row.names=F,quote=T,sep="\t")
    write.table(fullA,paste(outputname,".prm",sep=""),col.names=T,row.names=F,quote=T)
    cat(paste('....... Output file exported ',outputname,'.grm',' ..........\n',sep=''))
    cat(paste('....... Output file exported ',outputname,'.prm',' ..........\n',sep=''))
  }
  
  if(plots==TRUE){
    plotname <-paste("Genomic_relationship_",outputname,".png",sep="")
    dpi=300;width.cm<-15;height.cm<-13;width.in<-width.cm/2.54;height.in<-height.cm/2.54
    png(file=plotname,width=width.in*dpi,height=height.in*dpi,pointsize=8,units="px",res=dpi)
    par(mfrow=c(1,2))
    par(mar=c(5,4,5,2))
    br<-length(IID)*0.5; if(br>100){br=100}
    hist(diag(G),breaks=br,col=sample(2:6,1),xlab="diagonal element of G-matrix",
         main=paste("grm for ",outputname,"\n","based on ",ana_type,sep=""))
    hist(G[lower.tri(G)],breaks=br,col=sample(2:6,1),
         xlab="off-diagonal element of G-matrix",main="")
    dev.off()
  }
  return(G)
}
