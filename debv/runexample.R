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

source("dEBV.R")
debv <-dEBV(trait="Pheno_trait.txt",Pedigree="Pedigree_info.txt",genoIDs="IDs_genotyped.txt",
            dataformat="DIR",h2=0.30,p.varSNP=0.50,outname="GyrdEBV")



