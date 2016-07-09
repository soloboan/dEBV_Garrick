dEBV <- function(trait,Pedigree,genoIDs,dataformat,h2,p.varSNP,outname){
  cat('\n...........  Deregression procedure following Garrick et al. 2009  ...........\n')
  # Deregression of EBVs following Garrick et al. 2009
  # Part of the script were sourced from Badke et al. 2014
  
  #----------------------------------------------------------------------------------------------------------------------#
  # The deregression scripts requires the following parameters and information to run sucessfully
  
  #1. trait=""      +++++ The file containing the EBVs, columns = (ID, EBV, reliability)
  #2. Pedigree=""   +++++ The file containing the Pedigree information, columns = (ID, Sire, Dam)
  #3. genoIDs=""    +++++ The file containing the IDs of the genotyped individuals
  #4. dataformat="" +++++ Specify if the files are in the directory or it an R-object (options are 'DIR','R-object')
  #5. h2=""         +++++ heritabilty of the traits
  #6. p.varSNP=""   +++++ proportion of genetic variance accounted for by marker genotypes 
  #7. outname=""    +++++ Output filename 
  #----------------------------------------------------------------------------------------------------------------------#

  if (dataformat=="DIR"){
    Pedigree <- read.table(Pedigree,colClasses=c("character","character","character"))
    cat('\n...........  Pedigree file imported  ...........\n')
    trait <- read.table(trait,colClasses=c("character","numeric","numeric"))
    cat('...........  Phenotype file imported  ...........\n')
    genoIDs <- read.table(genoIDs,colClasses=c("character"))
    cat('...........  IDs of genotyped individuals imported  ...........\n')
  } else if (dataformat=="Robject"){
    Pedigree <- get(Pedigree)
    cat('\n...........  Pedigree file read  ...........\n')
    trait <- get(trait)
    cat('...........  Phenotype file read  ...........\n')
    genoIDs <- get(genoIDs)
    cat('...........  IDs of genotyped individuals read  ...........\n')
  }
  
  cat('\n...........  Preparing file for deregression  ...........\n')
  
  PedTraits <- merge(x=trait,y=Pedigree,by.x=1,by.y=1)
  PedTraitsgenoIDs <- merge(x=PedTraits,y=genoIDs,by.x=1,by.y=1)
  colnames(x=PedTraitsgenoIDs) <- c("ID","ID_EBV","ID_Rel","SireID","DamID") 
  Rel.sire <- merge(x=PedTraits,y=PedTraitsgenoIDs,by.x=1,by.y=4)
  Rel.sire <- Rel.sire[,-4:-5]  
  colnames(Rel.sire)[1:3] <- c('SireID','Sire_EBV','Sire_R2')
  Rel.dam <- merge(x=PedTraits,y=Rel.sire,by.x=1,by.y=7)
  Rel.dam <- Rel.dam[,-4:-5]  
  colnames(Rel.dam)[1:3] <- c('DamID','Dam_EBV','Dam_R2')
  data <- Rel.dam[,c('ID','ID_EBV','ID_Rel','SireID','Sire_EBV','Sire_R2','DamID','Dam_EBV','Dam_R2')]
  dEBV <- data[,-c(1,4,7)]
  dEBV$h2 <- h2
  dEBV$p.varSNP <- 1-p.varSNP
  
  Debv_Garrick <- function (ebv_mat){
    lambda = (1 - ebv_mat[7])/ebv_mat[7]
    PA = (ebv_mat[3] + ebv_mat[5])/2
    rPA = (ebv_mat[4] + ebv_mat[6])/4
    alpha = 1/(0.5 - rPA)
    delta = (0.5 - rPA)/(1 - ebv_mat[2])
    ZpZPA = lambda*(0.5*alpha - 4) + 0.5*lambda*sqrt(alpha^2 +16/delta)
    ZpZi = delta*ZpZPA + 2*lambda*(2*delta - 1)
    LHS = rbind(cbind(ZpZPA + 4*lambda, -2*lambda),cbind(-2*lambda, ZpZi + 2*lambda))
    #L1 = solve(LHS)
    RHS = LHS %*% c(PA, ebv_mat[1])
    drgi = RHS[2]/ZpZi
    rdrg = 1 - lambda/(ZpZi + lambda)
    we = (1 - ebv_mat[7])/((ebv_mat[8] + (1 - rdrg)/rdrg)*ebv_mat[7])
    ret = c(ebv_mat[1],ebv_mat[2],round(drgi,3),round(rdrg,3),round(we,3))
    names(ret) <- c("EBV", "r2EBV", "dEBV","r2dEBV","weight")
    return(ret)
  }
  
  cat('...........  File preparation finished  ...........\n\n')
  cat('...........  Deregression process started  ...........\n')
  cat('...........  Deregression script is a modified version of Badke et al. 2014  ...........\n')
  
  Deregress <- t(apply(X=dEBV,MARGIN=1,FUN=Debv_Garrick))
  Deregress <- cbind.data.frame(IID=data[,1],Deregress)
  write.table(Deregress,paste(outname,".debv",sep=""),col.names=T,row.names=F,quote=F,sep='\t')
  
  cat('...........  Deregression finished    ...........\n\n')
  cat(paste('...........  Output files exported as ',outname,'.debv ...........',sep=''),'\n')
  
  return(Deregress)
}

